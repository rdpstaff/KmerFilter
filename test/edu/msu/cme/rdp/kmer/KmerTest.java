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

import java.util.Random;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class KmerTest {

    @Test
    public void testKmer() {

        Random rand = new Random();

        for (int length = 1; length < 65; length++) {
            for (int iteration = 0; iteration < 10; iteration++) {
                StringBuilder strKmer = new StringBuilder();
                for (int index = 0; index < length; index++) {
                    strKmer.append(Kmer.intToChar[rand.nextInt(4)]);
                }

                Kmer kmer = new Kmer(strKmer.toString().toCharArray());
                assertEquals(strKmer.toString(), kmer.toString());

                StringBuilder rightShift = new StringBuilder(strKmer.toString());

                char shift = Kmer.intToChar[rand.nextInt(4)];
                rightShift.insert(0, shift);
                rightShift.deleteCharAt(rightShift.length() - 1);
                assertEquals("Shift Right Original kmer: " + strKmer + ", length: " + length + ", shift=" + shift,rightShift.toString(), kmer.shiftRight(shift).toString());

                StringBuilder leftShift = new StringBuilder(strKmer.toString());
                shift = Kmer.intToChar[rand.nextInt(4)];

                leftShift.append(shift);
                leftShift.deleteCharAt(0);
		assertEquals("Shift Left Original kmer: " + strKmer + ", length: " + length + ", shift=" + shift, leftShift.toString(), kmer.shiftLeft(shift).toString());
            }
        }
    }
}
