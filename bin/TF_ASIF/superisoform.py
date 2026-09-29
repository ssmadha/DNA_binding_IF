import warnings


class Segment:
    """
    One non-overlapping piece of a gene's coding region, as used to build
    a segment-based superisoform.

    Segments come from cutting every coding exon (CDS block) of every
    transcript at every CDS boundary found across the gene's transcripts,
    so two transcripts either share a segment exactly or not at all. A
    segment is identified by its genomic span and its phase: the same
    span read in a different frame encodes different residues, so it is
    a different segment.
    """

    def __init__(self, index, start, end, phase):
        """
        Constructor

        Parameters
        ----------
        index: int
            Position of this segment in the superisoform (5' to 3' in
            the gene's orientation).
        start: int
            Genomic start (1-based, inclusive).
        end: int
            Genomic end (1-based, inclusive).
        phase: int
            Bases to skip from this segment's 5' end (in transcript
            orientation) to reach the first codon starting inside it;
            same convention as the GTF "frame" column.
        """
        self.index = index
        self.start = start
        self.end = end
        self.phase = phase
        self.transcripts = set()
        self.nt_seq = None

    @property
    def key(self):
        return self.start, self.end, self.phase

    def __len__(self):
        return self.end - self.start + 1

    def __repr__(self):
        return "Segment %d %d-%d phase %d in %d transcripts" % (
            self.index, self.start, self.end, self.phase, len(self.transcripts))


class Superisoform:
    """
    Segment-based superisoform: the union of a gene's coding segments
    across all its transcripts, plus each transcript's path through them.
    """

    def __init__(self, gene_id, strand, transcript_blocks, cds_sequences=None):
        """
        Constructor

        Parameters
        ----------
        gene_id: str
            Ensembl Gene ID.
        strand: int
            1 or -1.
        transcript_blocks: dict[str, list[tuple]]
            For each Ensembl Transcript ID, its CDS blocks as (genomic
            start, genomic end, GTF frame) tuples in transcript (5' to
            3') order, excluding the stop codon.
        cds_sequences: dict[str, Seq], optional
            Ensembl CDS FASTA sequences keyed by transcript ID, used to
            fill in each segment's nucleotide sequence.
        """
        self.gene_id = gene_id
        self.strand = strand
        self.transcript_blocks = transcript_blocks
        # Leading Ns Ensembl's CDS FASTA adds to an incomplete (5'
        # truncated) CDS, per the GTF frame of its first block.
        self.padding = {enst_id: (3 - int(blocks[0][2])) % 3 if blocks[0][2].isdigit() else 0
                        for enst_id, blocks in transcript_blocks.items()}

        breakpoints = set()
        for blocks in transcript_blocks.values():
            for start, end, _ in blocks:
                breakpoints.add(start)
                breakpoints.add(end + 1)
        breakpoints = sorted(breakpoints)

        # (enst_id, [(piece start, piece end, phase, offset into CDS)]) in
        # transcript order
        pieces_by_transcript = {}
        for enst_id, blocks in transcript_blocks.items():
            pad = self.padding[enst_id]
            offset = 0
            pieces = []
            for start, end, _ in blocks:
                cuts = [bp for bp in breakpoints if start < bp <= end] + [end + 1]
                spans = []
                piece_start = start
                for cut in cuts:
                    spans.append((piece_start, cut - 1))
                    piece_start = cut
                if strand == -1:
                    spans = spans[::-1]
                for piece_start, piece_end in spans:
                    phase = (3 - (pad + offset) % 3) % 3
                    pieces.append((piece_start, piece_end, phase, offset))
                    offset += piece_end - piece_start + 1
            pieces_by_transcript[enst_id] = pieces

        keys = {(start, end, phase) for pieces in pieces_by_transcript.values()
                for start, end, phase, _ in pieces}
        ordered_keys = sorted(keys, key=lambda key: (key[0] * strand, key[2]))
        self.segments = [Segment(i, *key) for i, key in enumerate(ordered_keys)]
        segment_index = {segment.key: segment.index for segment in self.segments}

        self.transcript_segments = {}
        self._transcript_offsets = {}
        self._path_positions = {}
        for enst_id, pieces in pieces_by_transcript.items():
            indices = [segment_index[(start, end, phase)] for start, end, phase, _ in pieces]
            self.transcript_segments[enst_id] = indices
            self._transcript_offsets[enst_id] = [offset for _, _, _, offset in pieces]
            for index in indices:
                self.segments[index].transcripts.add(enst_id)

        if cds_sequences is not None:
            self._fill_sequences(cds_sequences)

    def _fill_sequences(self, cds_sequences):
        """
        Set each segment's nucleotide sequence from the CDS FASTA
        sequence of a transcript containing it, warning if two
        transcripts disagree on it.

        Parameters
        ----------
        cds_sequences: dict[str, Seq]
        """
        for enst_id, indices in self.transcript_segments.items():
            cds_seq = cds_sequences.get(enst_id)
            if cds_seq is None:
                continue
            pad = self.padding[enst_id]
            for index, offset in zip(indices, self._transcript_offsets[enst_id]):
                segment = self.segments[index]
                nt_seq = str(cds_seq[pad + offset:pad + offset + len(segment)])
                if segment.nt_seq is None:
                    segment.nt_seq = nt_seq
                elif segment.nt_seq != nt_seq:
                    warnings.warn("%s: segment %d sequence differs between transcripts (%s)"
                                  % (self.gene_id, index, enst_id))

    def transcript_nt_sequence(self, enst_id):
        """
        Rebuild a transcript's coding sequence (without padding or stop
        codon) by concatenating its segments' sequences.

        Parameters
        ----------
        enst_id: str

        Returns
        -------
        str
        """
        return "".join(self.segments[index].nt_seq for index in self.transcript_segments[enst_id])

    def residue_codons(self, enst_id):
        """
        Map each residue of a transcript's protein to the exact coding
        bases of its codon, as (segment index, offset within segment)
        pairs. Padding bases of an incomplete-5' CDS have no segment and
        are left out, so the first residue of such a CDS has fewer than
        three bases.

        Parameters
        ----------
        enst_id: str

        Returns
        -------
        list[tuple[tuple[int, int], ...]]
            Codon bases for each residue, in protein order. Covers every
            complete codon, so a complete CDS's stop codon (not part of
            the GTF CDS blocks) is not included.
        """
        per_base = [None] * self.padding[enst_id]
        for index in self.transcript_segments[enst_id]:
            per_base += [(index, offset) for offset in range(len(self.segments[index]))]
        return [tuple(base for base in per_base[codon_start:codon_start + 3] if base is not None)
                for codon_start in range(0, len(per_base) - 2, 3)]

    def residue_segments(self, enst_id):
        """
        Map each residue of a transcript's protein to the segment(s)
        its codon lies in (see residue_codons). A codon split across a
        segment boundary maps to two (or, for a 1-2 nt segment, three)
        segments.

        Parameters
        ----------
        enst_id: str

        Returns
        -------
        list[tuple[int, ...]]
            Segment indices for each residue, in protein order.
        """
        return [tuple(dict.fromkeys(index for index, _ in codon)) for codon in self.residue_codons(enst_id)]

    def contains_codon(self, enst_id, codon):
        """
        Whether a transcript encodes exactly this codon: it contains
        every segment the codon's bases lie in and, for a codon split
        across segments, those segments are consecutive in the
        transcript (otherwise it joins the codon's bases to a different
        exon and so encodes a different residue).

        Parameters
        ----------
        enst_id: str
        codon: tuple[tuple[int, int], ...]
            Codon bases, as returned by residue_codons.

        Returns
        -------
        bool
        """
        if enst_id not in self._path_positions:
            self._path_positions[enst_id] = {index: position for position, index
                                             in enumerate(self.transcript_segments[enst_id])}
        positions = self._path_positions[enst_id]
        indices = list(dict.fromkeys(index for index, _ in codon))
        if any(index not in positions for index in indices):
            return False
        return all(positions[b] == positions[a] + 1 for a, b in zip(indices, indices[1:]))

    def __repr__(self):
        return "Superisoform of %s: %d segments from %d transcripts" % (
            self.gene_id, len(self.segments), len(self.transcript_segments))