
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
package edu.msu.cme.rdp.kmer;

import java.util.Arrays;
import java.util.Random;

/**
 *
 * @author fishjord
 */
public class Kmer {

    public static final int a = 0, c = 1, g = 2, t = 3;
    public static final char[] intToChar = new char[]{'a', 'c', 'g', 't'};
    public static final byte[] validateLookup = new byte[127];
    public static final byte[] complementLookup = new byte[127];

    static {
        validateLookup['A'] = a;
        validateLookup['C'] = c;
        validateLookup['G'] = g;
        validateLookup['T'] = t;
        validateLookup['a'] = a;
        validateLookup['c'] = c;
        validateLookup['g'] = g;
        validateLookup['t'] = t;
        validateLookup['U'] = t;
        validateLookup['u'] = t;

        Arrays.fill(complementLookup, (byte) -1);

        complementLookup['A'] = t;  // A
        complementLookup['T'] = a;  // T
        complementLookup['U'] = a;  // U
        complementLookup['G'] = c;  // G
        complementLookup['C'] = g;  // C

        complementLookup[a] = t;  // a
        complementLookup[t] = a;  // t
        complementLookup[g] = c;  // g
        complementLookup[c] = g;  // c
    }
    private long kmer1;
    private long kmer2;
    private final int k;
    private final int l;
    private final int lastFill;
    private final long lastMask;

    /**
     *
     * @param charKmer
     */
    public Kmer(char[] charKmer) {
        this(charKmer.length);
        byte b;

        /*
         * So this is mildly complicated We want to store the kmer in
         * little-endian byte order so that the shifts make logical sense (right
         * = >> and left = <<), but that means storing the kmer kind of
         * backwards in the array...
         */

        for (int index = 0; index < charKmer.length; index++) {
            b = validateLookup[charKmer[index]];
            if (b == -1) {
                throw new IllegalArgumentException("Kmer contains one or more invalid bases (" + new String(charKmer) + ")");
            }

            if (index < 32) {
                kmer1 = (kmer1 << 2) | (b & 0x3);
            } else {
                kmer2 = (kmer2 << 2) | (b & 0x3);
            }
        }
    }

    public Kmer(long kmer1, int k) {
        this(kmer1, 0, k);
        if(k > 32) {
            throw new IllegalArgumentException("this constructor can only accept k sizes <= 32");
        }
    }

    public Kmer(long kmer1, long kmer2, int k) {
        this(k);
        this.kmer1 = kmer1;
        this.kmer2 = kmer2;
    }

    private Kmer(int k) {
        if (k > 64) {
            throw new IllegalArgumentException("k-mer lenght must be <= 64");
        }

        this.k = k;
        l = (int) Math.ceil(k / 32.0);
        lastFill = 32 - (32 * l - k);
        long tmp = 0;
        for (int index = 0; index < lastFill; index++) {
            tmp = (tmp << 2) | 3;
        }
        lastMask = tmp;
    }

    private Kmer(Kmer k) {
        this.k = k.k;
        this.l = k.l;
        this.kmer1 = k.kmer1;
        this.kmer2 = k.kmer2;

        this.lastFill = k.lastFill;
        this.lastMask = k.lastMask;
    }

    public int length() {
        return k;
    }

    public int packedLength() {
        return l;
    }

    public long getPart(int index) {
        if (index == 0) {
            return kmer1;
        }
        return kmer2;
    }

    public Kmer shiftRight(char c) {
        return shiftRight(validateLookup[c]);
    }

    public Kmer shiftLeft(char c) {
        return shiftLeft(validateLookup[c]);
    }

    public Kmer shiftRight(byte b) {
        Kmer ret = new Kmer(this);
        byte overflow = b;
        byte tmp;

        if (l == 2) {
            tmp = (byte) (kmer1 & 0x3);
            ret.kmer1 = (kmer1 >>> 2) | ((long) overflow << 62);
            overflow = tmp;

            ret.kmer2 = (kmer2 >>> 2) | ((long) overflow << ((ret.lastFill - 1) * 2));
            ret.kmer2 &= lastMask;
        } else {
            ret.kmer1 = (kmer1 >>> 2) | ((long) overflow << ((ret.lastFill - 1) * 2));
            ret.kmer1 &= lastMask;
        }

        return ret;
    }

    public Kmer shiftLeft(byte b) {
        Kmer ret = new Kmer(this);
        byte overflow = b;
        byte tmp;

        if (l == 2) {
            tmp = (byte) ((kmer2 >>> ((lastFill - 1) * 2)) & 0x3);
            ret.kmer2 = (kmer2 << 2) | ((long) overflow & 0x3);
            ret.kmer2 &= lastMask;
            overflow = tmp;

            ret.kmer1 = (kmer1 << 2) | (overflow);
        } else {
            ret.kmer1 = (kmer1 << 2) | ((long) overflow & 0x3);
            ret.kmer1 &= lastMask;
        }

        return ret;
    }

    @Override
    public String toString() {
        StringBuilder buf = new StringBuilder();
        for (int seg = l; seg > 0; seg--) {
            int index = 32;
            if (seg == l) {
                index = lastFill;
            }
            long thisk = (seg == 2)? kmer2 : kmer1;
            for (; index > 0; index--) {
                buf.append(intToChar[(int) (thisk & 0x3)]);
                thisk = thisk >> 2;
            }
        }

        return buf.reverse().toString();
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Kmer other = (Kmer) obj;
        if (this.kmer1 != other.kmer1 || this.kmer2 != other.kmer2) {
            return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 59 * hash + (int) (this.kmer1 ^ (this.kmer1 >>> 32));
        hash = 59 * hash + (int) (this.kmer2 ^ (this.kmer2 >>> 32));
        return hash;
    }

    public static Kmer randomKmer(int l) {
        Random rand = new Random();
        char[] ckmer = new char[l];
        for(int index = 0;index < ckmer.length;index++) {
            ckmer[index] = Kmer.intToChar[rand.nextInt(4)];
        }

        return new Kmer(ckmer);
    }
}