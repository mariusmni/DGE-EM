/*
 *  Copyright 2011 marius.
 * 
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 * 
 *       http://www.apache.org/licenses/LICENSE-2.0
 * 
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *  under the License.
 */

package edu.uconn.engr.dna.dgeem;

import java.util.Map;

/**
 *
 * @author marius
 */
public class WeightedTagClassJ {
	public String[] isoforms;
	public final int[] pos;
	public final double[] w;
	public int[] iso;
	public int multiplicity;

	public WeightedTagClassJ(String[] isoforms, int[] pos, double[] w,
			int multiplicity) {
		this.isoforms = isoforms;
		this.pos = pos;
		this.w = w;
		this.multiplicity = multiplicity;
//		assert(isoforms.length == pos.length && pos.length == w.length,
//		   "All arrays must have the same length; found " + isoforms.length +
//		   " " + pos.length + " " + w.length);
	}

	public int size() {
		return pos.length;
	}

	public void incrementMultiplicity() {
		multiplicity += 1;
	}

	// only called before running EM
	public void indexData(Map<String, Integer> index) {
		iso = new int[isoforms.length];
		for (int i = 0; i < iso.length; ++i) {
			iso[i] = index.get(isoforms[i]);
		}
 	}

	public String[] isoforms() {
		return isoforms;
	}
	
	@Override
	public boolean equals(Object o) {
		if (!(o instanceof WeightedTagClassJ)) {
			return false;
		}
		WeightedTagClassJ t = (WeightedTagClassJ)o;
        if (size() != t.size()) {
            return false;
		}

        for (int i = 0; i < isoforms.length; ++i) {
            if (pos[i] != t.pos[i] ||
				Double.doubleToLongBits(w[i]) != Double.doubleToLongBits(t.w[i]) ||
				!isoforms[i].equals(t.isoforms[i])) {
                return false;
			}
		}
		return true;
	}

	@Override
	public int hashCode() {
		int h = 0;
        for (int i = 0; i < isoforms.length; ++i) {
			h = 31 * h + pos[i];
            long bits = Double.doubleToLongBits(w[i]);
			h = 31 * h + (int)(bits ^ (bits >>> 32));
			String element = isoforms[i];
			h = 31 * h + (element == null ? 0 : element.hashCode());
		}
		return h;
	}
}





