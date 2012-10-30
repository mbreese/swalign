#!/usr/bin/env python
'''
Simple Smith-Waterman aligner
'''
import sys
import StringIO


class ScoringMatrix(object):
    '''
    Read scoring matrix from a file or string

    Matrix should be space-delimited in a format like:

      A C G T
    A 1 0 0 0
    C 0 1 0 0
    G 0 0 1 0
    T 0 0 0 1

    Rows and Columns must be in the same order

    '''
    def __init__(self, filename=None, text=None):
        assert filename or text

        if filename:
            fs = open(filename)
        else:
            fs = StringIO.StringIO(text)

        self.scores = []
        self.bases = None

        for line in fs:
            if line[0] == '#':
                continue

            if not self.bases:
                self.bases = line.split()
                self.base_count = len(self.bases)
            else:
                cols = line.split()
                self.scores.extend([int(x) for x in cols[1:]])

        fs.close()

    def score(self, one, two):
        one_idx = 0
        two_idx = 0
        for i, b in enumerate(self.bases):
            if b == one:
                one_idx = i
            if b == two:
                two_idx = i

        return self.scores[(one_idx * self.base_count) + two_idx]


class NucleotideScoringMatrix(object):
    def __init__(self, match=1, mismatch=-1):
        self.match = match
        self.mismatch = mismatch

    def score(self, one, two):
        if one == two:
            return self.match
        return self.mismatch


class Matrix(object):
    def __init__(self, rows, cols, init=None):
        self.rows = rows
        self.cols = cols
        self.values = [init, ] * rows * cols

    def get(self, row, col):
        return self.values[(row * self.cols) + col]

    def set(self, row, col, val):
        self.values[(row * self.cols) + col] = val


class LocalAlignment(object):
    def __init__(self, scoring_matrix, gap_penalty=-1, gap_extension_penalty=-1, gap_extension_decay=0.0, prefer_gap_runs=True, verbose=False):
        self.scoring_matrix = scoring_matrix
        self.gap_penalty = gap_penalty
        self.gap_extension_penalty = gap_extension_penalty
        self.gap_extension_decay = gap_extension_decay
        self.verbose = verbose
        self.prefer_gap_runs = prefer_gap_runs

    def align(self, ref, query, ref_name='', query_name='', rc=False):
        ref = ref.upper()
        query = query.upper()

        matrix = Matrix(len(query) + 1, len(ref) + 1, (0, ' ', 0))

        max_val = 0
        max_row = 0
        max_col = 0

        # calculate matrix
        for row in xrange(1, matrix.rows):
            for col in xrange(1, matrix.cols):
                mm_val = matrix.get(row - 1, col - 1)[0] + self.scoring_matrix.score(query[row - 1], ref[col - 1])

                ins_run = 0
                del_run = 0

                if matrix.get(row - 1, col)[1] == 'i':
                    ins_run = matrix.get(row - 1, col)[2]
                    if matrix.get(row - 1, col)[0] == 0:
                        # no penalty to start the alignment
                        ins_val = 0
                    else:
                        if not self.gap_extension_decay:
                            ins_val = matrix.get(row - 1, col)[0] + self.gap_extension_penalty
                        else:
                            ins_val = matrix.get(row - 1, col)[0] + min(0, self.gap_extension_penalty + ins_run * self.gap_extension_decay)
                else:
                    ins_val = matrix.get(row - 1, col)[0] + self.gap_penalty

                if matrix.get(row, col - 1)[1] == 'd':
                    del_run = matrix.get(row, col - 1)[2]
                    if matrix.get(row, col - 1)[0] == 0:
                        # no penalty to start the alignment
                        del_val = 0
                    else:
                        if not self.gap_extension_decay:
                            del_val = matrix.get(row, col - 1)[0] + self.gap_extension_penalty
                        else:
                            del_val = matrix.get(row, col - 1)[0] + min(0, self.gap_extension_penalty + del_run * self.gap_extension_decay)

                else:
                    del_val = matrix.get(row, col - 1)[0] + self.gap_penalty

                cell_val = max(mm_val, del_val, ins_val, 0)

                if not self.prefer_gap_runs:
                    ins_run = 0
                    del_run = 0

                if del_run and cell_val == del_val:
                    val = (cell_val, 'd', del_run + 1)
                elif ins_run and cell_val == ins_val:
                    val = (cell_val, 'i', ins_run + 1)
                elif cell_val == mm_val:
                    val = (cell_val, 'm', 0)
                elif cell_val == del_val:
                    val = (cell_val, 'd', 1)
                elif cell_val == ins_val:
                    val = (cell_val, 'i', 1)
                else:
                    val = (0, 'x', 0)

                if val[0] >= max_val:
                    max_val = val[0]
                    max_row = row
                    max_col = col

                matrix.set(row, col, val)

        # backtrack
        row = max_row
        col = max_col
        val = max_val
        op = 'm'

        aln = []

        path = []
        while val > 0:
            path.append((row, col))
            aln.append(op)

            op = matrix.get(row, col)[1]

            if op == 'm':
                row -= 1
                col -= 1
            elif op == 'i':
                row -= 1
            elif op == 'd':
                col -= 1
            else:
                break

            val, op, runlen = matrix.get(row, col)
        aln.reverse()

        if self.verbose:
            self.dump_matrix(ref, query, matrix, path)
            print aln
            print (max_row, max_col), max_val

        cigar = _reduce_cigar(aln)
        # if not self.regap:
        # else:
        #     cigar = []
        #     last_del = 0
        #     refpos = 0
        #     qpos = 0

        #     changed = False
        #     for size, op in _reduce_cigar(aln):
        #         if op == 'M':
        #             if last_del > size and cigar[-1][1] == 'M':
        #                 gap_seq = ref[col + refpos - last_del:col + refpos]
        #                 fragment = query[row + qpos: row + qpos + size]

        #                 adjustment = 0

        #                 for a, b in zip(fragment, gap_seq):
        #                     if a == b:
        #                         adjustment += 1
        #                     else:
        #                         break

        #                 if adjustment:
        #                     if self.verbose:
        #                         print "gap_seq :", gap_seq
        #                         print "fragment:", fragment
        #                         print "adjustment:", adjustment
        #                     size -= adjustment
        #                     cigar[-1] = (cigar[-1][0] + adjustment, 'M')
        #                     changed = True

        #             if last_del:
        #                 cigar.append((last_del, 'D'))
        #                 last_del = 0

        #             refpos += size
        #             qpos += size

        #             cigar.append((size, op))
        #         elif op == 'D':
        #             refpos += size
        #             last_del = size
        #         elif op == 'I':
        #             qpos += size
        #             if last_del:
        #                 cigar.append((last_del, 'D'))
        #                 last_del = 0

        #             cigar.append((size, op))

        #     if last_del:
        #         cigar.append((last_del, 'D'))

        #     if self.verbose and changed:
        #         print "re-aligning deletes"
        #         print "old cigar: ", _cigar_str(_reduce_cigar(aln))
        #         print "new cigar: ", _cigar_str(cigar)
        #         Alignment(query, ref, row, col, _reduce_cigar(aln), max_val, ref_name, query_name, rc).dump()

        return Alignment(query, ref, row, col, cigar, max_val, ref_name, query_name, rc)

    def dump_matrix(self, ref, query, matrix, path, show_row=-1, show_col=-1):
        sys.stdout.write('      -      ')
        sys.stdout.write('       '.join(ref))
        sys.stdout.write('\n')
        for row in xrange(matrix.rows):
            if row == 0:
                sys.stdout.write('-')
            else:
                sys.stdout.write(query[row - 1])

            for col in xrange(matrix.cols):
                if show_row == row and show_col == col:
                    sys.stdout.write('       *')
                else:
                    sys.stdout.write(' %5s%s%s' % (matrix.get(row, col)[0], matrix.get(row, col)[1], '$' if (row, col) in path else ' '))
            sys.stdout.write('\n')


def _reduce_cigar(operations):
    count = 1
    last = None
    ret = []
    for op in operations:
        if last and op == last:
            count += 1
        elif last:
            ret.append((count, last.upper()))
            count = 1
        last = op

    if last:
        ret.append((count, last.upper()))
    return ret


def _cigar_str(cigar):
    out = ''
    for num, op in cigar:
        out += '%s%s' % (num, op)
    return out


class Alignment(object):
    def __init__(self, query, ref, q_pos, r_pos, cigar, score, ref_name='', query_name='', rc=False):
        self.query = query
        self.ref = ref
        self.q_pos = q_pos
        self.r_pos = r_pos
        self.cigar = cigar
        self.score = score
        self.r_name = ref_name
        self.q_name = query_name
        self.rc = rc

        q_len = 0
        r_len = 0

        self.matches = 0
        self.mismatches = 0

        i = self.r_pos
        j = self.q_pos

        for count, op in self.cigar:
            if op == 'M':
                q_len += count
                r_len += count
                for k in xrange(count):
                    if self.query[j] == self.ref[i]:
                        self.matches += 1
                    else:
                        self.mismatches += 1
                    i += 1
                    j += 1

            elif op == 'I':
                q_len += count
                j += count
                self.mismatches += count
            elif op == 'D':
                r_len += count
                i += count
                self.mismatches += count

        self.q_end = q_pos + q_len
        self.r_end = r_pos + r_len
        if self.mismatches + self.matches > 0:
            self.identity = float(self.matches) / (self.mismatches + self.matches)
        else:
            self.identity = 0

    @property
    def extended_cigar_str(self):
        qpos = 0
        rpos = 0
        ext_cigar_str = ''
        working = []
        for count, op in self.cigar:
            if op == 'M':
                for k in xrange(count):
                    if self.query[self.q_pos + qpos + k] == self.ref[self.r_pos + rpos + k]:
                        ext_cigar_str += 'M'
                    else:
                        ext_cigar_str += 'X'
                qpos += count
                rpos += count

            elif op == 'I':
                qpos += count
                ext_cigar_str += 'I' * count
            elif op == 'D':
                rpos += count
                ext_cigar_str += 'D' * count

            working = _reduce_cigar(ext_cigar_str)

        out = ''
        for num, op in working:
            out += '%s%s' % (num, op)
        return out

    @property
    def cigar_str(self):
        return _cigar_str(self.cigar)

    def dump(self, out=sys.stdout):
        i = self.r_pos
        j = self.q_pos

        if not self.rc:
            q = 'Query: %4s ' % self.q_pos
        else:
            q = 'Query: %4s ' % (len(self.query) - self.q_pos)
        r = 'Ref  : %4s ' % self.r_pos
        m = '            '

        q_len = 0
        for count, op in self.cigar:
            if op == 'M':
                q_len += count
                for k in xrange(count):
                    q += self.query[j]
                    r += self.ref[i]
                    if self.query[j] == self.ref[i]:
                        m += '|'
                    else:
                        m += '.'

                    i += 1
                    j += 1
            elif op == 'D':
                for k in xrange(count):
                    q += '-'
                    r += self.ref[i]
                    m += ' '
                    i += 1
            elif op == 'I':
                q_len += count
                for k in xrange(count):
                    q += self.query[j]
                    r += '-'
                    m += ' '
                    j += 1

            elif op == 'N':
                q += '-//-'
                r += '-//-'
                m += '    '

        if self.q_name:
            out.write('Query: %s%s (%s nt)\n' % (self.q_name, ' (reverse-compliment)' if self.rc else '', len(self.query)))
        if self.r_name:
            out.write('Ref  : %s (%s nt)\n\n' % (self.r_name, len(self.ref)))

        if not self.rc:
            out.write("%s %s\n" % (q, self.q_end))
        else:
            out.write("%s %s\n" % (q, len(self.query) - self.q_pos - q_len + 1))

        out.write("%s\n" % m)
        out.write("%s %s\n" % (r, self.r_end))
        out.write("Score: %s\n" % self.score)
        out.write("Matches: %s (%.1f%%)\n" % (self.matches, self.identity * 100))
        out.write("Mismatches: %s\n" % (self.mismatches,))
        out.write("CIGAR: %s\n" % self.cigar_str)


def fasta_gen(fname):
    def gen():
        seq = ''
        name = ''

        if fname == '-':
            f = sys.stdin
            name = 'stdin'
        else:
            f = open(fname)

        for line in f:
            if line[0] == '>':
                if name and seq:
                    yield (name, seq)

                name = line[1:].strip().split(' ')[0]
                seq = ''
            else:
                seq += line.strip()

        if name and seq:
            yield (name, seq)

        if fname != '-':
            f.close()
    return gen


def seq_gen(name, seq):
    def gen():
        yield (name, seq)

    return gen

__revcomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
__cache = {}


def revcomp(seq):
    if seq in __cache:
        return __cache[seq]

    ret = []
    for s in seq.upper()[::-1]:
        ret.append(__revcomp[s])

    __cache[seq] = ''.join(ret)
    return __cache[seq]


#     sw.align('ACACACTA','AGCACACA').dump()
#     aln=sw.align("AAGGGGAGGACGATGCGGATGTTC","AGGGAGGACGATGCGG")
#     aln.dump()
