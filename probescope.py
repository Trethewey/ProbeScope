#!/usr/bin/env python3
"""
ProbeScope
==========
Generates per-gene coverage track plots from BAM files overlaid on
panel target BED annotations with MANE-Select exon models.

Usage examples:
  # Basic: one BAM, one target BED, one probe BED
  python probescope.py \
    --bams sample1.bam \
    --probes probe_design.bed \
    --targets panel_targets.bed \
    --output results/

  # Two BAMs, two target BEDs, specific genes
  python probescope.py \
    --bams sample1.bam sample2.bam \
    --bam-labels "Sample A" "Sample B" \
    --probes CAM_Lymphoma.bed \
    --targets gaea_panel.bed athena_panel.bed \
    --target-labels "GAEA Panel" "Athena Panel" \
    --genes TP53 BRCA1 MYC \
    --output results/

  # All genes in target BED (no gap filtering)
  python probescope.py \
    --bams sample.bam \
    --probes probes.bed \
    --targets targets.bed \
    --all-genes \
    --output results/

Requirements:
  pip install bamnostic pandas numpy matplotlib

No reference genome database needed - MANE-Select exon coordinates
are fetched from the UCSC API and cached locally as JSON.

BED file formats supported:
  - 4-column: chrom  start  end  transcript
  - 6-column: chrom  start  end  transcript  gene  exon
  - Chromosome names with or without 'chr' prefix (auto-detected)

"""

import argparse
import json
import sys
import urllib.request
import urllib.error
from pathlib import Path

import bamnostic as bs
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter


# ═══════════════════════════════════════════════════════════════════════════════
# BED PARSING
# ═══════════════════════════════════════════════════════════════════════════════

def detect_chr_prefix(filepath):
    """Check whether a BED file uses 'chr' prefix on chromosome names."""
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts and parts[0]:
                return parts[0].startswith('chr')
    return False


def parse_bed(filepath, transcript_to_gene=None):
    """
    Parse a BED file. Auto-detects format:
      4-col: chrom start end transcript
      6-col: chrom start end transcript gene exon
    Normalises chromosome names to bare format (no 'chr' prefix).
    """
    rows = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('track'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue

            chrom = parts[0].replace('chr', '')
            start = int(parts[1])
            end = int(parts[2])

            transcript = parts[3].split(';')[0] if len(parts) >= 4 else ''
            gene = parts[4] if len(parts) >= 5 else ''
            exon = parts[5] if len(parts) >= 6 else ''

            # If no gene name in BED, try mapping from transcript
            if not gene and transcript and transcript_to_gene:
                gene = transcript_to_gene.get(transcript, transcript)
            elif not gene:
                gene = transcript

            rows.append({
                'chrom': chrom, 'start': start, 'end': end,
                'transcript': transcript, 'gene': gene, 'exon': exon,
            })

    df = pd.DataFrame(rows)
    if len(df) > 0:
        df['size'] = df['end'] - df['start']
    return df


def build_transcript_to_gene(bed_dfs):
    """Build a transcript-to-gene mapping from any BED that has gene names."""
    t2g = {}
    for df in bed_dfs:
        if 'transcript' in df.columns and 'gene' in df.columns:
            for _, row in df.iterrows():
                if row['transcript'] and row['gene'] and row['gene'] != row['transcript']:
                    t2g[row['transcript']] = row['gene']
    return t2g


# ═══════════════════════════════════════════════════════════════════════════════
# GAP COMPUTATION
# ═══════════════════════════════════════════════════════════════════════════════

def subtract_intervals(target_start, target_end, covers):
    """Return portions of [target_start, target_end) not covered by any interval in covers."""
    uncovered = [(target_start, target_end)]
    for cs, ce in sorted(covers):
        new = []
        for us, ue in uncovered:
            if ce <= us or cs >= ue:
                new.append((us, ue))
            else:
                if cs > us:
                    new.append((us, cs))
                if ce < ue:
                    new.append((ce, ue))
        uncovered = new
    return uncovered


def compute_gaps(target_df, probe_by_chrom):
    """Find regions in target_df not covered by probe regions."""
    gap_rows = []
    for _, g in target_df.iterrows():
        covers = probe_by_chrom.get(g['chrom'], [])
        relevant = [(s, e) for s, e in covers if s < g['end'] and e > g['start']]
        for us, ue in subtract_intervals(g['start'], g['end'], relevant):
            gap_rows.append({
                'chrom': g['chrom'], 'start': us, 'end': ue,
                'gene': g['gene'], 'transcript': g['transcript'],
                'gap_size': ue - us,
            })
    return pd.DataFrame(gap_rows) if gap_rows else pd.DataFrame()


# ═══════════════════════════════════════════════════════════════════════════════
# MANE-SELECT EXON DATA
# ═══════════════════════════════════════════════════════════════════════════════

def load_mane_cache(cache_path):
    """Load cached MANE exon data from JSON file."""
    if cache_path.exists():
        with open(cache_path) as f:
            return json.load(f)
    return {}


def save_mane_cache(cache_path, mane_data):
    """Save MANE exon data to JSON cache."""
    with open(cache_path, 'w') as f:
        json.dump(mane_data, f, indent=2)


def fetch_mane_exons(gene, transcript, chrom_bare, region_start, region_end):
    """Fetch MANE-Select exon coordinates from UCSC API for a single gene."""
    chrom = f"chr{chrom_bare}"
    padded_start = max(0, region_start - 50000)
    padded_end = region_end + 50000

    url = (f"https://api.genome.ucsc.edu/getData/track?genome=hg38"
           f"&track=ncbiRefSeq&chrom={chrom}&start={padded_start}&end={padded_end}")

    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=20) as resp:
            data = json.loads(resp.read().decode())

        track_key = 'ncbiRefSeq' if 'ncbiRefSeq' in data else list(data.keys())[0]
        if isinstance(data.get(track_key), list):
            # Try exact transcript match first, then prefix match
            transcript_base = transcript.rsplit('.', 1)[0] if transcript else ''
            for item in data[track_key]:
                name = item.get('name', '')
                if name.startswith(transcript_base):
                    exon_starts = [int(x) for x in item['exonStarts'].rstrip(',').split(',') if x]
                    exon_ends = [int(x) for x in item['exonEnds'].rstrip(',').split(',') if x]
                    return {
                        'chrom': chrom,
                        'strand': item.get('strand', '+'),
                        'txStart': int(item['txStart']),
                        'txEnd': int(item['txEnd']),
                        'exonStarts': exon_starts,
                        'exonEnds': exon_ends,
                    }
    except Exception as e:
        print(f"    [!!] Failed to fetch MANE data for {gene}: {e}", file=sys.stderr)

    return None


def ensure_mane_data(genes_info, mane_cache, cache_path):
    """
    Ensure MANE exon data exists for all requested genes.
    Fetches missing data from UCSC API and updates the cache.

    genes_info: list of dicts with keys: gene, transcript, chrom, start, end
    """
    missing = [g for g in genes_info if g['gene'] not in mane_cache]
    if not missing:
        return mane_cache

    print(f"  Fetching MANE-Select exon data for {len(missing)} genes from UCSC API...")
    for g in missing:
        result = fetch_mane_exons(g['gene'], g['transcript'], g['chrom'],
                                  g['start'], g['end'])
        if result:
            mane_cache[g['gene']] = result
            n_exons = len(result['exonStarts'])
            print(f"    [OK] {g['gene']}: {n_exons} exons ({result['strand']})")
        else:
            print(f"    [--] {g['gene']}: not found")

    save_mane_cache(cache_path, mane_cache)
    return mane_cache


# ═══════════════════════════════════════════════════════════════════════════════
# BAM DEPTH
# ═══════════════════════════════════════════════════════════════════════════════

def get_depth_array(bam_handle, chrom, start, end):
    """
    Build a depth array for a genomic region from BAM.
    Does a single index lookup and iterates reads once.
    Returns (positions, depths) arrays downsampled for plotting.
    """
    size = end - start
    if size <= 0:
        return np.zeros(0), np.zeros(0)

    depth_array = np.zeros(size, dtype=np.int32)
    try:
        for read in bam_handle.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                continue
            rs = read.reference_start
            re = read.reference_end if read.reference_end else rs + read.query_length
            cs = max(rs, start) - start
            ce = min(re, end) - start
            if ce > cs:
                depth_array[cs:ce] += 1
    except Exception:
        pass

    # Downsample for plotting (max 3000 points)
    if size > 3000:
        n_bins = 3000
        positions = np.linspace(start, end - 1, n_bins, dtype=int)
        bin_size = max(1, size // n_bins)
        depths = np.array([
            depth_array[max(0, p - start):min(p - start + bin_size, size)].mean()
            for p in positions
        ])
    else:
        positions = np.arange(start, end)
        depths = depth_array.astype(float)

    return positions, depths


# ═══════════════════════════════════════════════════════════════════════════════
# PLOTTING
# ═══════════════════════════════════════════════════════════════════════════════

# Colour palettes for multiple targets/BAMs
TARGET_COLOURS = [
    ('#7E57C2', '#4527A0'),  # purple
    ('#26A69A', '#00695C'),  # teal
    ('#FF7043', '#D84315'),  # deep orange
    ('#5C6BC0', '#283593'),  # indigo
    ('#66BB6A', '#2E7D32'),  # green
]

BAM_COLOURS = [
    '#1565C0',  # blue
    '#C62828',  # red
    '#2E7D32',  # green
    '#F57F17',  # amber
    '#6A1B9A',  # purple
]

GAP_COLOURS = [
    ('#F44336', '#B71C1C'),  # red
    ('#FF6F00', '#E65100'),  # orange
    ('#AB47BC', '#6A1B9A'),  # purple
    ('#42A5F5', '#1565C0'),  # blue
    ('#26A69A', '#00695C'),  # teal
]


def plot_gene(gene, target_dfs, target_labels, target_gaps,
              probe_df, bam_handles, bam_labels, mane_data, ax_cov, ax_tracks):
    """Plot coverage + annotation tracks for one gene."""

    mane = mane_data.get(gene)

    # Determine x-axis range from all data sources
    all_starts, all_ends = [], []
    for tdf in target_dfs:
        g = tdf[tdf['gene'] == gene]
        if len(g) > 0:
            all_starts.extend(g['start'].tolist())
            all_ends.extend(g['end'].tolist())

    if len(probe_df) > 0:
        chrom_bare = None
        for tdf in target_dfs:
            g = tdf[tdf['gene'] == gene]
            if len(g) > 0:
                chrom_bare = g['chrom'].iloc[0]
                break
        if chrom_bare:
            p = probe_df[(probe_df['chrom'] == chrom_bare)]
            # Narrow to gene locus
            if all_starts and all_ends:
                p = p[(p['start'] < max(all_ends) + 5000) & (p['end'] > min(all_starts) - 5000)]
            if len(p) > 0:
                all_starts.extend(p['start'].tolist())
                all_ends.extend(p['end'].tolist())

    if mane:
        all_starts.append(mane['txStart'])
        all_ends.append(mane['txEnd'])

    if not all_starts:
        return

    global_min = min(all_starts)
    global_max = max(all_ends)
    span = global_max - global_min
    pad = span * 0.03

    # Get chromosome
    chrom_bare = None
    for tdf in target_dfs:
        g = tdf[tdf['gene'] == gene]
        if len(g) > 0:
            chrom_bare = g['chrom'].iloc[0]
            break

    if not chrom_bare:
        return

    # ── Coverage plot ──
    fetch_start = max(0, global_min - int(pad))
    fetch_end = global_max + int(pad)
    print(f"    Fetching depth for {gene} ({chrom_bare}:{fetch_start:,}-{fetch_end:,})...",
          end=" ", flush=True)

    for bi, (bam_h, bam_label) in enumerate(zip(bam_handles, bam_labels)):
        colour = BAM_COLOURS[bi % len(BAM_COLOURS)]
        pos, depth = get_depth_array(bam_h, chrom_bare, fetch_start, fetch_end)
        if len(pos) > 0:
            ax_cov.fill_between(pos, depth, alpha=0.35, color=colour, linewidth=0)
            ax_cov.plot(pos, depth, color=colour, linewidth=0.5, alpha=0.7, label=bam_label)

    print("done", flush=True)

    # Shade gap regions
    for gi, gaps in enumerate(target_gaps):
        gene_gaps = gaps[gaps['gene'] == gene] if len(gaps) > 0 else pd.DataFrame()
        if len(gene_gaps) > 0:
            gc = GAP_COLOURS[gi % len(GAP_COLOURS)][0]
            for _, g in gene_gaps.iterrows():
                ax_cov.axvspan(g['start'], g['end'], alpha=0.12, color=gc, zorder=0)

    # Coverage axis formatting
    ax_cov.set_xlim(global_min - pad, global_max + pad)
    y_max = ax_cov.get_ylim()[1]
    if y_max < 1:
        y_max = 100
    ax_cov.set_ylim(0, y_max * 1.15)
    ax_cov.set_ylabel('Read Depth', fontsize=8)
    ax_cov.axhline(y=30, color='#FF9800', linewidth=0.8, linestyle='--', alpha=0.6)
    ax_cov.text(global_max + pad * 0.3, 30, '30x', fontsize=6, color='#FF9800', va='center')
    ax_cov.legend(fontsize=7, loc='upper right', framealpha=0.8)
    ax_cov.grid(axis='y', alpha=0.15, linestyle='-')
    ax_cov.yaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{x:,.0f}'))
    for spine in ['top', 'right']:
        ax_cov.spines[spine].set_visible(False)

    # Title
    n_exons = len(mane['exonStarts']) if mane else '?'
    strand_str = f' ({mane["strand"]})' if mane else ''
    probe_bp = probe_df[
        (probe_df['chrom'] == chrom_bare) &
        (probe_df['start'] < global_max + 5000) &
        (probe_df['end'] > global_min - 5000)
    ]['size'].sum() if len(probe_df) > 0 else 0
    gap_strs = []
    for gi, (gaps, label) in enumerate(zip(target_gaps, target_labels)):
        gene_gaps = gaps[gaps['gene'] == gene] if len(gaps) > 0 else pd.DataFrame()
        gbp = gene_gaps['gap_size'].sum() if len(gene_gaps) > 0 else 0
        if gbp > 0:
            gap_strs.append(f'{label} gap: {gbp:,}')
    gap_info = ' | '.join(gap_strs) if gap_strs else 'no gaps'
    ax_cov.set_title(
        f'{gene}  --  {n_exons} exons{strand_str} | Probes: {probe_bp:,} bp | {gap_info}',
        fontsize=10, fontweight='bold', loc='left'
    )

    # ── Annotation tracks ──
    n_targets = len(target_dfs)
    # Tracks: MANE, targets (1-N), probes, gaps
    n_tracks = 1 + n_targets + 1 + 1  # MANE + targets + probes + gaps
    th = 0.22
    sg = 0.10

    y_positions = {}
    y_idx = n_tracks - 1
    y_positions['MANE-Select'] = y_idx * (th + sg)
    y_idx -= 1
    for label in target_labels:
        y_positions[label] = y_idx * (th + sg)
        y_idx -= 1
    y_positions['Probes'] = y_idx * (th + sg)
    y_idx -= 1
    y_positions['Gaps'] = y_idx * (th + sg)

    # MANE exons
    if mane:
        y = y_positions['MANE-Select']
        mid = y + th / 2
        ax_tracks.plot([mane['txStart'], mane['txEnd']], [mid, mid],
                       color='#333', linewidth=1.2, zorder=1)
        arr = '>' if mane['strand'] == '+' else '<'
        n_arr = max(1, int(span / (span * 0.04)))
        for apos in np.linspace(mane['txStart'] + span*0.02, mane['txEnd'] - span*0.02, n_arr):
            ax_tracks.text(apos, mid, arr, ha='center', va='center',
                          fontsize=5, color='#888', fontweight='bold', zorder=0)
        for ie, (es, ee) in enumerate(zip(mane['exonStarts'], mane['exonEnds'])):
            ax_tracks.add_patch(plt.Rectangle(
                (es, y), ee-es, th,
                facecolor='#37474F', edgecolor='#263238', linewidth=0.6, alpha=0.9, zorder=2))
            if (ee-es) > span * 0.012:
                ax_tracks.text((es+ee)/2, y+th+0.02, f'E{ie+1}',
                              ha='center', va='bottom', fontsize=5, color='#37474F', fontweight='bold')

    # Target tracks
    for ti, (tdf, label) in enumerate(zip(target_dfs, target_labels)):
        y = y_positions[label]
        fc, ec = TARGET_COLOURS[ti % len(TARGET_COLOURS)]
        gene_t = tdf[tdf['gene'] == gene]
        for _, r in gene_t.iterrows():
            ax_tracks.add_patch(plt.Rectangle(
                (r['start'], y), r['end']-r['start'], th,
                facecolor=fc, edgecolor=ec, linewidth=0.4, alpha=0.75))

    # Probe track
    y_probe = y_positions['Probes']
    probe_gene = probe_df[
        (probe_df['chrom'] == chrom_bare) &
        (probe_df['start'] < global_max + 5000) &
        (probe_df['end'] > global_min - 5000)
    ]
    for _, r in probe_gene.iterrows():
        ax_tracks.add_patch(plt.Rectangle(
            (r['start'], y_probe), r['end']-r['start'], th,
            facecolor='#2196F3', edgecolor='#1565C0', linewidth=0.4, alpha=0.85))

    # Gap tracks (stacked within one row)
    y_gap = y_positions['Gaps']
    n_gap_tracks = max(1, len(target_gaps))
    sub_h = th / n_gap_tracks
    for gi, (gaps, label) in enumerate(zip(target_gaps, target_labels)):
        gene_gaps = gaps[gaps['gene'] == gene] if len(gaps) > 0 else pd.DataFrame()
        fc, ec = GAP_COLOURS[gi % len(GAP_COLOURS)]
        for _, g in gene_gaps.iterrows():
            ax_tracks.add_patch(plt.Rectangle(
                (g['start'], y_gap + gi * sub_h), g['end']-g['start'], sub_h * 0.9,
                facecolor=fc, edgecolor=ec, linewidth=0.6, alpha=0.85))

    # Track axis formatting
    ax_tracks.set_xlim(global_min - pad, global_max + pad)
    ax_tracks.set_ylim(-0.08, n_tracks * (th + sg) + 0.08)
    yticks = [y + th/2 for y in y_positions.values()]
    ax_tracks.set_yticks(yticks)
    ax_tracks.set_yticklabels(list(y_positions.keys()), fontsize=7, fontweight='bold')
    ax_tracks.set_xlabel(f'Genomic Position (chr{chrom_bare})', fontsize=8)
    ax_tracks.xaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{x:,.0f}'))
    ax_tracks.grid(axis='x', alpha=0.1, linestyle='-')
    ax_tracks.set_axisbelow(True)
    for spine in ['top', 'right']:
        ax_tracks.spines[spine].set_visible(False)


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        prog='probescope',
        description='ProbeScope: Generate per-gene coverage track plots from BAM files '
                    'overlaid on panel target BED annotations with MANE-Select exon models.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # One BAM, one target BED
  python probescope.py --bams sample.bam --probes probes.bed --targets panel.bed

  # Two BAMs, two target BEDs, custom labels
  python probescope.py \\
    --bams sample1.bam sample2.bam \\
    --bam-labels "Sample A" "Sample B" \\
    --probes CAM_probes.bed \\
    --targets gaea.bed athena.bed \\
    --target-labels "GAEA" "Athena" \\
    --output results/

  # Specific genes only
  python probescope.py --bams s.bam --probes p.bed --targets t.bed --genes TP53 MYC BCL2

  # All genes (not just those with gaps)
  python probescope.py --bams s.bam --probes p.bed --targets t.bed --all-genes
        """
    )
    parser.add_argument('--bams', nargs='+', required=True,
                        help='One or more BAM files (must have .bai index)')
    parser.add_argument('--bam-labels', nargs='+', default=None,
                        help='Labels for BAM files (default: filenames)')
    parser.add_argument('--probes', required=True,
                        help='Probe design BED file (the regions your assay targets)')
    parser.add_argument('--targets', nargs='+', required=True,
                        help='One or more panel target BED files to compare against probes')
    parser.add_argument('--target-labels', nargs='+', default=None,
                        help='Labels for target BED files (default: filenames)')
    parser.add_argument('--genes', nargs='+', default=None,
                        help='Specific gene names to plot (default: all genes with gaps)')
    parser.add_argument('--all-genes', action='store_true',
                        help='Plot all genes, not just those with gaps')
    parser.add_argument('--output', '-o', default='coverage_output',
                        help='Output directory (default: coverage_output/)')
    parser.add_argument('--mane-cache', default=None,
                        help='Path to MANE exon cache JSON (default: <output>/mane_cache.json)')
    parser.add_argument('--per-page', type=int, default=4,
                        help='Number of genes per page (default: 4)')
    parser.add_argument('--dpi', type=int, default=150,
                        help='Output DPI (default: 150)')

    args = parser.parse_args()

    # ── Validate inputs ────────────────────────────────────────────────────
    for bam_path in args.bams:
        if not Path(bam_path).exists():
            sys.exit(f"Error: BAM file not found: {bam_path}")
        bai = Path(bam_path).with_suffix('.bam.bai')
        if not bai.exists():
            # Try alternative index naming
            bai2 = Path(str(bam_path) + '.bai')
            if not bai2.exists():
                sys.exit(f"Error: BAM index not found: {bai} or {bai2}")

    if not Path(args.probes).exists():
        sys.exit(f"Error: Probe BED not found: {args.probes}")

    for t in args.targets:
        if not Path(t).exists():
            sys.exit(f"Error: Target BED not found: {t}")

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Labels ─────────────────────────────────────────────────────────────
    bam_labels = args.bam_labels or [Path(b).stem.split('_')[0] for b in args.bams]
    if len(bam_labels) != len(args.bams):
        sys.exit(f"Error: {len(bam_labels)} BAM labels provided but {len(args.bams)} BAM files")

    target_labels = args.target_labels or [Path(t).stem[:30] for t in args.targets]
    if len(target_labels) != len(args.targets):
        sys.exit(f"Error: {len(target_labels)} target labels provided but {len(args.targets)} target files")

    # ── Load BED files ─────────────────────────────────────────────────────
    print("Loading BED files...")

    # Parse all targets first to build transcript-to-gene map
    target_dfs_raw = [parse_bed(t) for t in args.targets]
    t2g = build_transcript_to_gene(target_dfs_raw)

    # Re-parse targets that lacked gene names, now with the mapping
    target_dfs = []
    for t, df_raw in zip(args.targets, target_dfs_raw):
        needs_remap = (df_raw['gene'] == df_raw['transcript']).all()
        if needs_remap and t2g:
            target_dfs.append(parse_bed(t, t2g))
        else:
            target_dfs.append(df_raw)

    probe_df = parse_bed(args.probes, t2g)

    for label, df in zip(target_labels, target_dfs):
        print(f"  {label}: {len(df)} regions, {df['size'].sum():,} bp, {df['gene'].nunique()} genes")
    print(f"  Probes: {len(probe_df)} regions, {probe_df['size'].sum():,} bp")

    # ── Compute gaps ───────────────────────────────────────────────────────
    print("\nComputing gaps...")
    probe_by_chrom = {}
    for _, r in probe_df.iterrows():
        probe_by_chrom.setdefault(r['chrom'], []).append((r['start'], r['end']))

    target_gaps = []
    for label, tdf in zip(target_labels, target_dfs):
        gaps = compute_gaps(tdf, probe_by_chrom)
        target_gaps.append(gaps)
        gap_bp = gaps['gap_size'].sum() if len(gaps) > 0 else 0
        gap_genes = gaps['gene'].nunique() if len(gaps) > 0 else 0
        print(f"  {label}: {len(gaps)} gap regions, {gap_bp:,} bp across {gap_genes} genes")

    # ── Determine genes to plot ────────────────────────────────────────────
    if args.genes:
        genes_to_plot = args.genes
    elif args.all_genes:
        all_genes = set()
        for tdf in target_dfs:
            all_genes |= set(tdf['gene'].unique())
        genes_to_plot = sorted(all_genes)
    else:
        # Only genes with gaps
        gap_genes = set()
        for gaps in target_gaps:
            if len(gaps) > 0:
                gap_genes |= set(gaps['gene'].unique())
        if not gap_genes:
            print("\nNo gaps found - all target regions are covered by probes.")
            print("Use --all-genes to plot all genes anyway.")
            sys.exit(0)
        # Sort by total gap size
        gene_gap_sizes = {}
        for gaps in target_gaps:
            if len(gaps) > 0:
                for gene in gaps['gene'].unique():
                    gene_gap_sizes[gene] = gene_gap_sizes.get(gene, 0) + gaps[gaps['gene'] == gene]['gap_size'].sum()
        genes_to_plot = sorted(gap_genes, key=lambda g: -gene_gap_sizes.get(g, 0))

    print(f"\n  Will plot {len(genes_to_plot)} genes")

    # ── MANE-Select data ───────────────────────────────────────────────────
    cache_path = Path(args.mane_cache) if args.mane_cache else output_dir / "mane_cache.json"
    mane_data = load_mane_cache(cache_path)

    genes_info = []
    for gene in genes_to_plot:
        for tdf in target_dfs:
            g = tdf[tdf['gene'] == gene]
            if len(g) > 0:
                genes_info.append({
                    'gene': gene,
                    'transcript': g['transcript'].iloc[0],
                    'chrom': g['chrom'].iloc[0],
                    'start': g['start'].min(),
                    'end': g['end'].max(),
                })
                break

    mane_data = ensure_mane_data(genes_info, mane_data, cache_path)

    # ── Open BAMs ──────────────────────────────────────────────────────────
    print("\nOpening BAM files...")
    bam_handles = []
    for bam_path in args.bams:
        bam_handles.append(bs.AlignmentFile(str(bam_path), 'rb'))
    print(f"  {len(bam_handles)} BAMs opened")

    # ── Generate plots ─────────────────────────────────────────────────────
    plt.rcParams.update({
        'font.size': 9, 'axes.titlesize': 11, 'axes.labelsize': 9,
        'figure.dpi': args.dpi, 'figure.facecolor': 'white',
    })

    page_size = args.per_page
    n_pages = (len(genes_to_plot) + page_size - 1) // page_size

    print(f"\nGenerating plots ({n_pages} pages, {page_size} genes per page)...")

    for page_idx in range(n_pages):
        page_genes = genes_to_plot[page_idx * page_size : (page_idx + 1) * page_size]
        n_g = len(page_genes)

        fig, axes = plt.subplots(n_g * 2, 1, figsize=(22, 6.5 * n_g),
                                 gridspec_kw={'height_ratios': [2, 1] * n_g})
        if n_g == 1:
            axes = list(axes)

        for i, gene in enumerate(page_genes):
            print(f"  [{page_idx+1}/{n_pages}] {gene}...")
            plot_gene(gene, target_dfs, target_labels, target_gaps,
                     probe_df, bam_handles, bam_labels, mane_data,
                     axes[i*2], axes[i*2+1])

        # Global legend
        legend_elements = []
        for bi, label in enumerate(bam_labels):
            legend_elements.append(mpatches.Patch(
                facecolor=BAM_COLOURS[bi % len(BAM_COLOURS)], alpha=0.4, label=f'Depth: {label}'))
        legend_elements.append(mpatches.Patch(
            facecolor='#37474F', edgecolor='#263238', label='MANE-Select Exons'))
        for ti, label in enumerate(target_labels):
            fc, ec = TARGET_COLOURS[ti % len(TARGET_COLOURS)]
            legend_elements.append(mpatches.Patch(facecolor=fc, edgecolor=ec, label=label))
        legend_elements.append(mpatches.Patch(
            facecolor='#2196F3', edgecolor='#1565C0', label='Probes'))
        for gi, label in enumerate(target_labels):
            fc, ec = GAP_COLOURS[gi % len(GAP_COLOURS)]
            legend_elements.append(mpatches.Patch(facecolor=fc, edgecolor=ec, label=f'{label} Gap'))

        fig.legend(handles=legend_elements, loc='upper center',
                   ncol=min(len(legend_elements), 8), fontsize=8,
                   bbox_to_anchor=(0.5, 1.003), frameon=True, fancybox=True)

        page_label = f" (page {page_idx+1}/{n_pages})" if n_pages > 1 else ""
        fig.suptitle(
            f'ProbeScope{page_label}',
            fontsize=13, fontweight='bold', y=1.02
        )
        plt.tight_layout()
        suffix = f"_p{page_idx+1}" if n_pages > 1 else ""
        outpath = output_dir / f"probescope_output{suffix}.png"
        fig.savefig(outpath, dpi=args.dpi, bbox_inches='tight')
        print(f"  [OK] Saved: {outpath}")
        plt.close()

    # ── Save summary CSV ───────────────────────────────────────────────────
    summary_rows = []
    for gene in genes_to_plot:
        row = {'gene': gene}
        for label, tdf, gaps in zip(target_labels, target_dfs, target_gaps):
            gt = tdf[tdf['gene'] == gene]
            gg = gaps[gaps['gene'] == gene] if len(gaps) > 0 else pd.DataFrame()
            row[f'{label}_total_bp'] = gt['size'].sum() if len(gt) > 0 else 0
            row[f'{label}_gap_bp'] = gg['gap_size'].sum() if len(gg) > 0 else 0
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    summary_path = output_dir / "gap_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"\n  [OK] Summary: {summary_path}")

    print(f"\nDone. All outputs in: {output_dir}")


if __name__ == '__main__':
    main()
