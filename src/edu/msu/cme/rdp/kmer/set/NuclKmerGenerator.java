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

import java.util.Iterator;

/**
 *
 * @author fishjord
 */
public class NuclKmerGenerator
        implements KmerGenerator {

    private char[] bases;
    private int index;
    private int position;
    private long value;
    private Long next;
    private final long mask;
    private final int k;
    private boolean modelOnly;

    public NuclKmerGenerator(String seq, int k) {
        this(seq, k, false);
    }

    public NuclKmerGenerator(String seq, int k, boolean modelOnly) {
        if (k > 31) {
            throw new IllegalArgumentException("K-mer size cannot be larger than 31");
        }

        /*if (seq.length() < k) {
            throw new IllegalArgumentException("Sequence length is less than the kmer length");
        }*/

        bases = seq.toCharArray();
        mask = 0xffffffffffffffffL ^ (0xffffffffffffffffL << (k * 2));
        this.k = k;
        this.modelOnly = modelOnly;

        value = 0;
        index = 0;
        position = 1;
        next = nextInternal(0);
    }

    private void push(char c) {
        byte add;
        switch (c) {
            case 'a':
            case 'A':
                add = 0;
                break;
            case 'c':
            case 'C':
                add = 1;
                break;
            case 'g':
            case 'G':
                add = 2;
                break;
            case 't':
            case 'T':
            case 'u':
            case 'U':
                add = 3;
                break;
            default:
                throw new IllegalArgumentException("Only nucleotide bases excepted, not '" + c + "'");
        }
        value = value << 2;
        value |= add;
        value &= mask;
    }

    private Long nextInternal() {
        return nextInternal(k - 1);
    }

    private Long nextInternal(int klength) {
        while (index < bases.length) {
            char base = bases[index++];

            if (modelOnly && (Character.isLowerCase(base) || base == '-')) {
                if(base == '-') {
                    position++;
                }
                klength = 0;
            } else if(base == 'N' || base == 'n') {
                klength = 0;
            } else {
                if (base != '.') {
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
            ret.append(tochar((int) (v & 0x3)));
            v >>= 2;
        }

        return ret.reverse().toString();
    }

    private static final char tochar(int c) {
        switch (c) {
            case 0:
                return 'a';
            case 1:
                return 'c';
            case 2:
                return 'g';
            case 3:
                return 't';
            default:
                return '?';
        }
    }
}
