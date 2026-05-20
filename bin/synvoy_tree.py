"""Shared phylogeny utilities for SynVoy plots.

This module is intentionally dependency-free (no ete3, no numpy) so it can be
imported from any of the bin/ scripts without coupling to the iqtree+ete3
toolchain. ete3 trees can still be passed to `partition_clades` etc. via the
small `from_ete3` adapter.

Public surface:

    TreeNode                 — minimal node with parent links
    parse_newick_tree(text)  — newick → TreeNode root (None on empty)
    parse_newick_leaf_order(text)  — newick → leaf names left→right
    midpoint_root(root)      — re-root at midpoint of longest leaf-to-leaf path
    partition_clades(root, depth_frac=0.4)
                             — leaf name → clade_id; depth-threshold partition
    CLADE_PALETTE            — 12 distinct, color-blind-friendly hex colours
    species_from_leaf(name)  — `GOI_X|Sp_genus_fna_b0...` → `Sp_genus`
"""

import re
from collections import OrderedDict
from typing import Optional


# ─────────────────────────────── TreeNode ────────────────────────────────────

class TreeNode:
    """Lightweight phylogeny node with parent links."""
    __slots__ = ("name", "dist", "children", "parent")

    def __init__(self, name: str = "", dist: float = 0.0):
        self.name = name
        self.dist = float(dist or 0.0)
        self.children: list = []
        self.parent: Optional["TreeNode"] = None

    def is_leaf(self) -> bool:
        return not self.children

    def add(self, child: "TreeNode") -> None:
        self.children.append(child)
        child.parent = self

    def leaves(self):
        if self.is_leaf():
            yield self
            return
        for c in self.children:
            yield from c.leaves()

    def all_nodes(self):
        yield self
        for c in self.children:
            yield from c.all_nodes()


# ─────────────────────────────── parsing ─────────────────────────────────────

def parse_newick_tree(text: str) -> Optional[TreeNode]:
    """Parse a newick string into a TreeNode root. Returns None if empty."""
    if not text:
        return None
    s = text.strip().rstrip(";").strip()
    if not s:
        return None
    pos = [0]

    def _parse_node(parent=None):
        node = TreeNode()
        if parent is not None:
            parent.add(node)
        if pos[0] < len(s) and s[pos[0]] == "(":
            pos[0] += 1
            _parse_node(node)
            while pos[0] < len(s) and s[pos[0]] == ",":
                pos[0] += 1
                _parse_node(node)
            if pos[0] < len(s) and s[pos[0]] == ")":
                pos[0] += 1
        # name (or bootstrap label) until : , ) ;
        name_buf = []
        while pos[0] < len(s) and s[pos[0]] not in ":,();":
            name_buf.append(s[pos[0]])
            pos[0] += 1
        node.name = _strip_quotes("".join(name_buf))
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            num_buf = []
            while pos[0] < len(s) and s[pos[0]] not in ",();":
                num_buf.append(s[pos[0]])
                pos[0] += 1
            try:
                node.dist = float("".join(num_buf))
            except ValueError:
                node.dist = 0.0
        return node

    return _parse_node()


def parse_newick_leaf_order(text: str) -> list:
    """leaf names in left-to-right order. Empty list if parse fails."""
    root = parse_newick_tree(text)
    if root is None:
        return []
    return [n.name for n in root.leaves() if n.name]


def _strip_quotes(s: str) -> str:
    s = s.strip()
    if s.startswith("'") and s.endswith("'"):
        s = s[1:-1]
    return s


# ─────────────────────────────── midpoint rooting ───────────────────────────

def _patristic(a: TreeNode, b: TreeNode) -> float:
    """Distance between two leaves via their LCA. O(n) with parent links."""
    seen = {}
    n = a
    cum = 0.0
    while n is not None:
        seen[id(n)] = cum
        cum += n.dist
        n = n.parent
    n = b
    cum_b = 0.0
    while n is not None:
        if id(n) in seen:
            return seen[id(n)] + cum_b
        cum_b += n.dist
        n = n.parent
    return float("inf")


def midpoint_root(root: Optional[TreeNode]) -> Optional[TreeNode]:
    """Re-root the tree at the midpoint of the longest leaf-to-leaf path.

    Returns a NEW TreeNode (the input is left intact, since rerooting
    mutates parent links). If the tree has fewer than 2 leaves, returns
    the input unchanged.
    """
    if root is None:
        return None
    leaves = list(root.leaves())
    if len(leaves) < 2:
        return root

    def _farthest(start):
        best = (start, 0.0)
        for leaf in leaves:
            if leaf is start:
                continue
            d = _patristic(start, leaf)
            if d > best[1]:
                best = (leaf, d)
        return best

    a, _ = _farthest(leaves[0])
    b, total = _farthest(a)
    if total <= 0:
        return root
    half = total / 2.0

    # Path a → root and b → root.
    path_a = []
    n = a
    while n is not None:
        path_a.append(n)
        n = n.parent
    path_b = []
    n = b
    while n is not None:
        path_b.append(n)
        n = n.parent
    set_a = {id(x): i for i, x in enumerate(path_a)}
    lca = None
    i_lca = 0
    j_lca = 0
    for j, x in enumerate(path_b):
        if id(x) in set_a:
            lca = x
            i_lca = set_a[id(x)]
            j_lca = j
            break
    if lca is None:
        return root

    sequence = path_a[:i_lca + 1] + path_b[:j_lca][::-1]
    cum = 0.0
    for k in range(len(sequence) - 1):
        u = sequence[k]
        v = sequence[k + 1]
        if u.parent is v:
            edge = u.dist
        elif v.parent is u:
            edge = v.dist
        else:
            edge = abs(u.dist) + abs(v.dist)
        if cum + edge >= half:
            offset = half - cum
            child = u if u.parent is v else (v if v.parent is u else u)
            return _reroot_on_edge(root, child, offset, total_edge=edge)
        cum += edge
    return root


def _reroot_on_edge(orig_root: TreeNode, child: TreeNode,
                    offset_from_child_side: float, total_edge: float) -> TreeNode:
    """Build a new rooted tree where the root sits on the child-parent edge.

    `offset_from_child_side` units away from `child`. The child's branch is
    split into two: child→new_root (offset) and new_root→other_side (rest).
    """
    new_root = TreeNode(name="", dist=0.0)
    side_a = _clone_subtree(child)
    side_a.dist = offset_from_child_side
    new_root.add(side_a)
    other_side = _build_inverted_side(
        orig_root, child, remaining=total_edge - offset_from_child_side)
    new_root.add(other_side)
    return new_root


def _clone_subtree(node: TreeNode) -> TreeNode:
    copy = TreeNode(name=node.name, dist=node.dist)
    for c in node.children:
        copy.add(_clone_subtree(c))
    return copy


def _build_inverted_side(orig_root: TreeNode, broken_child: TreeNode,
                         remaining: float) -> TreeNode:
    """Walk broken_child.parent up to orig_root, inverting parent links so each
    former-parent becomes a node hanging off the new root branch.
    """
    chain = []
    n = broken_child.parent
    while n is not None:
        chain.append(n)
        n = n.parent
    if not chain:
        return TreeNode(name="", dist=remaining)
    new_inverted_root = TreeNode(name=chain[0].name, dist=remaining)
    for c in chain[0].children:
        if c is broken_child:
            continue
        new_inverted_root.add(_clone_subtree(c))
    current = new_inverted_root
    for next_anc, next_anc_in_chain in zip(chain, chain[1:]):
        anc_node = TreeNode(name=next_anc_in_chain.name, dist=next_anc.dist)
        for c in next_anc_in_chain.children:
            if c is next_anc:
                continue
            anc_node.add(_clone_subtree(c))
        current.add(anc_node)
        current = anc_node
    return new_inverted_root


# ─────────────────────────────── species collapse ───────────────────────────

def collapse_to_one_leaf_per_species(root: Optional[TreeNode],
                                     species_of) -> Optional[TreeNode]:
    """Return a copy of ``root`` with at most one leaf per species.

    For each species (as identified by ``species_of(leaf.name)``), keep the
    first leaf encountered in left→right traversal order; prune the rest.
    Internal nodes left without leaves are removed; nodes left with a single
    surviving child are contracted (the child's branch length absorbs the
    contracted parent's branch length, so cumulative depth from the root is
    preserved).

    Leaves whose ``species_of(name)`` returns None are dropped — they cannot
    be aligned to a per-species matrix row anyway. Returns None if pruning
    empties the tree.
    """
    if root is None:
        return None

    keep_ids = set()
    seen_species = set()
    for leaf in root.leaves():
        sp = species_of(leaf.name)
        if sp is None or sp in seen_species:
            continue
        seen_species.add(sp)
        keep_ids.add(id(leaf))

    def _clone(n):
        if n.is_leaf():
            if id(n) not in keep_ids:
                return None
            return TreeNode(name=n.name, dist=n.dist)
        new_children = []
        for c in n.children:
            cc = _clone(c)
            if cc is not None:
                new_children.append(cc)
        if not new_children:
            return None
        if len(new_children) == 1:
            child = new_children[0]
            child.dist = (child.dist or 0.0) + (n.dist or 0.0)
            return child
        copy = TreeNode(name=n.name, dist=n.dist)
        for cc in new_children:
            copy.add(cc)
        return copy

    return _clone(root)


# ─────────────────────────────── clade partition ─────────────────────────────

def partition_clades(
    root: Optional[TreeNode],
    target_k: Optional[int] = 4,
    depth_frac: float = 0.4,
) -> dict:
    """Split the (already midpoint-rooted) tree into clades.

    Two modes — picked by which kwarg you set:

    1. ``target_k`` (default 4) — *topology-driven iterative split*.
       Start with the whole tree as a single clade and repeatedly split
       the *largest* clade into its immediate children until ``K`` clades
       exist. Splits that would produce a singleton (one-leaf) child are
       deprioritised: if there's any multi-leaf-only split available we
       take that first. This matches what a human reader would draw
       around major lineages (e.g. ants vs. bumblebees+wasps vs. other
       bees), not "the one weird long-branch paralog vs. everyone else".

       Branch lengths are intentionally NOT used here. IQ-TREE outputs
       often have huge variation in leaf-branch length from divergent
       paralogs, and that variation dominates any length-based criterion
       — the longest branches are almost always outlier *leaves*, not
       internal splits between lineages.

    2. ``depth_frac`` (legacy) — *cumulative-depth threshold*.
       A node N is a clade root iff ``depth(parent) < threshold ≤ depth(N)``
       with ``threshold = depth_frac * max(depth(leaf))``. Kept for
       backward compatibility; not recommended.

    Returns ``{leaf_name: clade_id}`` with ``clade_id`` ≥ 0 assigned in
    left-to-right traversal order so adjacent leaves with the same id form
    visually contiguous groups in the rendered cladogram.
    """
    if root is None:
        return {}
    leaves = [n for n in root.all_nodes() if n.is_leaf()]
    if not leaves:
        return {}

    if target_k is not None:
        return _partition_top_down(root, leaves, int(target_k))
    return _partition_depth_threshold(root, leaves, float(depth_frac))


def _assign_clade_ids_lr(root: TreeNode, clade_root_per_leaf: dict) -> dict:
    """Helper: map clade-root → 0-based id in left-to-right leaf order."""
    seen = OrderedDict()
    out = {}
    for leaf in root.leaves():
        cr = clade_root_per_leaf.get(leaf.name)
        if cr is None:
            continue
        if cr not in seen:
            seen[cr] = len(seen)
        out[leaf.name] = seen[cr]
    return out


def _partition_top_down(root: TreeNode, leaves: list, k: int) -> dict:
    """Iteratively split the largest existing clade into its immediate
    children until we have ``k`` clades. Prefers splits whose children
    all have ≥ 2 leaves, so a long-branch outlier paralog doesn't get
    peeled off as its own one-leaf clade ahead of meaningful lineage
    boundaries.

    Midpoint rooting (the recommended pre-step) can produce *unary*
    internal nodes — single-child stubs left over from re-rooting on the
    middle of an edge. We follow those chains transparently so a stub
    never blocks a useful split.
    """
    if k <= 1:
        return {l.name: 0 for l in leaves}

    def n_leaves(n):
        return sum(1 for _ in n.leaves())

    def child_min_size(n):
        if not n.children:
            return 0
        return min(n_leaves(c) for c in n.children)

    def effective(n):
        """Descend any unary chain to the first real branching node (or leaf)."""
        while not n.is_leaf() and len(n.children) == 1:
            n = n.children[0]
        return n

    # Each "clade" is represented by a subtree-root node — but always its
    # *effective* root, i.e. the first node down a unary chain that has
    # ≥ 2 children or is a leaf. That way we can always ask "what are this
    # clade's two children?" without tripping on a single-child stub.
    #
    # Polytomies (n-ary internal nodes — common in NCBI taxonomy trees,
    # e.g. the Hymenoptera root has Apidae+Vespidae+Formicidae as three
    # siblings) need care: greedily expanding all children of the biggest
    # clade can blow past the K target in a single step. Instead we peel
    # children off the biggest clade one at a time, biggest-child first,
    # until either we've reached K or that clade has fully unfolded.
    #
    # When we peel off only some children of a polytomy, the un-peeled
    # children stay grouped under the original clade node. The
    # `clade_residual` map records which children belong to a "residual"
    # clade — needed because clade_node.leaves() returns *all* leaves below
    # it (including the peeled ones), and we want only the un-peeled.
    clades = [effective(root)]
    clade_residual = {}  # id(clade_node) → list of children to count as "this clade"
    while len(clades) < k:
        splittable = [c for c in clades if len(c.children) >= 2]
        if not splittable:
            break  # leaves only — nothing left to cut
        # First-class candidates: splits where every child has ≥ 2 leaves.
        non_singleton = [c for c in splittable if child_min_size(c) >= 2]
        pool = non_singleton if non_singleton else splittable
        biggest = max(pool, key=lambda c: (n_leaves(c), child_min_size(c)))

        # If we've previously peeled children off this node, only consider
        # its remaining residual children for this split.
        children = clade_residual.get(id(biggest)) or list(biggest.children)
        # Sort by size descending so the largest "peelable" group splits
        # off first — gives the most balanced cut at each step.
        children = sorted(children, key=n_leaves, reverse=True)

        # NB: compute `needed` BEFORE removing `biggest`. Net change per
        # iteration is +P (peel-and-keep-residual: remove biggest -1,
        # add residual +1, add P peeled = +P). So P = K - current_count.
        needed = k - len(clades)
        clades.remove(biggest)
        clade_residual.pop(id(biggest), None)

        if needed >= len(children) - 1:
            # Expand the whole polytomy: no residual, every child becomes
            # its own clade. Net change = len(children) - 1.
            clades.extend(effective(c) for c in children)
            continue

        # Peel off `needed` largest children; the rest stay grouped as a
        # single residual clade rooted at `biggest`.
        peeled = children[:needed]
        residual = children[needed:]
        for c in peeled:
            clades.append(effective(c))
        clade_residual[id(biggest)] = residual
        clades.append(biggest)

    # Build leaf → clade-root map. For residual clades, only count leaves
    # under the residual children, not the peeled ones.
    clade_root_per_leaf = {}
    for clade_node in clades:
        residual_children = clade_residual.get(id(clade_node))
        if residual_children is not None:
            for child in residual_children:
                for leaf in child.leaves():
                    clade_root_per_leaf[leaf.name] = id(clade_node)
        else:
            for leaf in clade_node.leaves():
                clade_root_per_leaf[leaf.name] = id(clade_node)
    return _assign_clade_ids_lr(root, clade_root_per_leaf)


def _partition_depth_threshold(root: TreeNode, leaves: list,
                               depth_frac: float) -> dict:
    """Legacy depth-fraction partition (kept for backward compatibility)."""
    depths = {}
    def _walk(node, d):
        depths[id(node)] = d
        for c in node.children:
            _walk(c, d + (c.dist or 0.0))
    _walk(root, 0.0)
    max_depth = max(depths[id(n)] for n in leaves)
    if max_depth <= 0:
        return {l.name: 0 for l in leaves}
    threshold = depth_frac * max_depth

    clade_root_per_leaf = {}
    for leaf in leaves:
        n = leaf
        chosen = leaf
        while n is not None:
            if depths[id(n)] >= threshold:
                chosen = n
            else:
                break
            n = n.parent
        clade_root_per_leaf[leaf.name] = id(chosen)

    return _assign_clade_ids_lr(root, clade_root_per_leaf)


# ─────────────────────────────── taxonomy tree ──────────────────────────────

def build_taxonomy_tree_newick(species_list, allow_network: bool = True):
    """Build a species tree from NCBI taxonomy (via ete3.NCBITaxa).

    Args:
        species_list: scientific names like 'Apis mellifera' or 'Apis_mellifera'.
            Underscores are tolerated and normalised to spaces for the lookup.
        allow_network: if False, only an already-downloaded taxdump is used
            (no first-run network fetch).

    Returns:
        Newick string with one leaf per resolvable species (name with
        underscores so it matches the matrix's species keys). Branch
        lengths are unit-1. Returns None if ete3 is missing or no species
        resolved — caller falls back to the GOI tree.

    The resulting tree captures phylogenetic relationships derived from
    classification alone (Apis_mellifera ↔ Apis_cerana ↔ Apis_florea
    share a genus node, etc.), which makes the matrix's tree column show
    who-is-related-to-who rather than per-locus sequence similarity.
    """
    try:
        from ete3 import NCBITaxa  # type: ignore
    except ImportError:
        return None

    try:
        if allow_network:
            ncbi = NCBITaxa()
        else:
            import os as _os
            default_db = _os.path.expanduser("~/.etetoolkit/taxa.sqlite")
            if not _os.path.exists(default_db):
                return None
            ncbi = NCBITaxa(dbfile=default_db)
    except Exception:
        return None

    # Normalise + dedupe input names while preserving caller order.
    clean = []
    seen = set()
    for sp in species_list:
        if not sp:
            continue
        norm = sp.replace("_", " ").strip()
        if norm and norm not in seen:
            seen.add(norm)
            clean.append(norm)
    if not clean:
        return None

    try:
        name2taxid = ncbi.get_name_translator(clean)
    except Exception:
        return None
    if not name2taxid:
        return None

    # Build taxid → species map. Prefer the first taxid ete3 returns.
    taxid_to_sp = {}
    for sp in clean:
        ids = name2taxid.get(sp)
        if ids:
            taxid_to_sp[ids[0]] = sp

    if len(taxid_to_sp) < 2:
        return None  # Need at least 2 leaves to form a tree.

    try:
        topo = ncbi.get_topology(list(taxid_to_sp.keys()),
                                 intermediate_nodes=False)
    except Exception:
        return None

    # Rename leaves from taxid → species (underscored), strip internal labels.
    for leaf in topo.iter_leaves():
        try:
            tid = int(leaf.name)
            leaf.name = taxid_to_sp.get(tid, leaf.name).replace(" ", "_")
        except (ValueError, KeyError):
            pass
    for nd in topo.traverse("postorder"):
        if not nd.is_leaf():
            nd.name = ""
        # NCBI topology has no real distances. Give every non-root node a
        # unit-1 branch so the matrix renderer can lay out the cladogram
        # horizontally; otherwise every internal node collapses to x=0.
        if nd.up is not None:
            nd.dist = 1.0

    # format=5 = newick with leaf names + branch lengths, no internal labels.
    return topo.write(format=5)


# Distinct-hue palette (Okabe-Ito + a few extras), enough for ~12 clades.
CLADE_PALETTE = [
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9",
    "#E69F00", "#882255", "#117733", "#999933", "#AA4499",
    "#44AA99", "#332288",
]


def color_for_clade(clade_id: int) -> str:
    return CLADE_PALETTE[clade_id % len(CLADE_PALETTE)]


# ─────────────────────────────── species helpers ─────────────────────────────

_SPECIES_TAIL_SUFFIXES = ("_exon_ann", "_fallback", "_flank_ann", "_full")


def species_from_leaf(leaf: str) -> Optional[str]:
    """`GOI_Melt|Apis_florea_fna_b0_l1_exon_ann` → `Apis_florea`.

    Also strips the ``_extraN`` suffix attached by the iterative-search
    tandem-copy renamer in ``iterative_search_runner.py`` (e.g.
    ``Apis_florea_fna_b0_l1_extra3`` → ``Apis_florea``). Without this,
    every renamed copy looks like its own species and downstream
    species-level collapses leak duplicate rows into the matrix
    cladogram.
    """
    if not leaf or "|" not in leaf:
        return None
    tail = leaf.split("|", 1)[1]
    for suf in _SPECIES_TAIL_SUFFIXES:
        if tail.endswith(suf):
            tail = tail[: -len(suf)]
    # Strip optional _extraN tag (tandem-copy disambiguator), then block/locus
    # IDs (_b0/_l1/_fl2 in any order), then the genome-extension stub.
    tail = re.sub(r"_extra\d+$", "", tail)
    tail = re.sub(r"_(b|fl|l)\d+(_(b|fl|l)\d+)*$", "", tail)
    tail = re.sub(r"_(fa|fna|fasta)$", "", tail)
    return tail
