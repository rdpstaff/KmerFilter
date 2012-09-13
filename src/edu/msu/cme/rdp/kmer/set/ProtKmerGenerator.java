/*
 * Copyright (C) 2012 Jordan Fish <fishjord at msu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package edu.msu.cme.rdp.kmer.set;

import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import java.util.Arrays;
import java.util.Iterator;

/**
 *
 * @author fishjord
 */
public class ProtKmerGenerator implements KmerGenerator {

    private char[] bases;
    private int index;
    private int position;
    private Long next = null;
    private long value;
    private final long mask;
    private final int k;
    private static final byte[] asciiMap = new byte[127];
    private boolean modelOnly = false;

    static {
        byte alpha = 0;
        Arrays.fill(asciiMap, (byte) -1);

        for (char c : SeqUtils.proteinAlphabet) {
            if (Character.isUpperCase(c)) {
                asciiMap[c] = asciiMap[Character.toLowerCase(c)] = alpha++;
            }
        }

        asciiMap['*'] = alpha++;

        if (alpha != 26) {
            throw new IllegalStateException("More than 25 amino acids...");
        }
    }

    public ProtKmerGenerator(String seq, int k) {
        this(seq, k, false);
    }

    public ProtKmerGenerator(String seq, int k, boolean modelOnly) {
        if (k > 10) {
            throw new IllegalArgumentException("K-mer size cannot be larger than 31");
        }

        if (seq.length() < k) {
            throw new IllegalArgumentException("Sequence length is less than the kmer length");
        }

        bases = seq.toCharArray();
        mask = 0xffffffffffffffffL ^ (0xffffffffffffffffL << (k * 5));
        this.k = k;
        this.modelOnly = modelOnly;

        value = 0;
        index = 0;
        position = 1;
        next = nextInternal(0);
    }

    private void push(char c) {
        if (asciiMap[c] == -1) {
            throw new IllegalArgumentException("Unknown prot base " + c);
        }
        value = value << 5;
        value |= asciiMap[c];
        value &= mask;
    }

    private Long nextInternal() {
        return nextInternal(k - 1);
    }

    private Long nextInternal(int klength) {
        while (index < bases.length) {
            char base = bases[index++];

            if (modelOnly && (Character.isLowerCase(base) || base == '-' || base == 'X' || base == 'x')) {
                if(base == '-' || base == 'X') {
                    position++;
                }
                klength = 0;
            } else {
                if (!modelOnly || (modelOnly && base != '.')) {
                    push(base);
                    position++;
                    klength++;
                }
            }
            if (klength == k) {
                return value;
            }
        }

        return null;
    }

    /**
     * Get the position of the FIRST CHARACTER in this kmer
     * in the sequence
     * @return
     */
    public int getPosition() {
        return position - k;
    }

    @Override
    public boolean hasNext() {
        if (next == null) {
            next = nextInternal();
        }

        return next != null;
    }

    @Override
    public Long next() {
        if (!hasNext()) {
            return null;
        }
        Long l = next;
        next = null;

        return l;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public String decodeLong(long v) {

        StringBuilder ret = new StringBuilder();
        for (int i = 0; i < k; i++) {
            ret.append(tochar((byte) (v & 0x1f)));
            v >>= 5;
        }

        return ret.reverse().toString();
    }

    private static final char tochar(byte c) {
        for (int index = 96; index < asciiMap.length; index++) {
            if (asciiMap[index] == c) {
                return (char) index;
            }
        }

        return '?';
    }

    public static void main(String[] args) {
        String seq = "tmamrqcalygkggigkstttqnlvaalaemgkkvmivgcdpkadstrlilhskaqgtvmemaasagsvedleledvlqigfggvkcvesggpepgvgcagrgvitainfleeegaysddldfvfydvlgdvvcg";
        ProtKmerGenerator kmer = new ProtKmerGenerator(seq, 5);
        int i = 0;
        while (kmer.hasNext()) {
            long v = kmer.next();
            String s = kmer.decodeLong(v);

            for (int index = 0; index < i; index++) {
                System.out.print(" ");
            }
            System.out.println(s);
            i++;
        }
    }
}
