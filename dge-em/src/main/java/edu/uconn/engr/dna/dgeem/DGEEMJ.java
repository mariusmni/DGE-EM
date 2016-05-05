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

import edu.cornell.lassp.houle.RngPack.Ranlux;
import edu.uconn.engr.dna.format.Isoforms;
import edu.uconn.engr.dna.util.Utils;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author marius
 */
public class DGEEMJ {
	private static final int maxSteps = 100;
	private static final double EPS = 0.00000000000000000001;
	private static final Object lock = new Object();
	public static int noConverge = 0;

	private Isoforms isoforms;
	private Map<String, Integer> nrCleavageSites;
	private Ranlux random = new Ranlux(4, 123);

	public DGEEMJ(Isoforms isoforms, Map<String, Integer> nrCleavageSites) {
		this.isoforms = isoforms;
		this.nrCleavageSites = nrCleavageSites;
	}


	private void println(String s) {
		System.out.println(s);
	}

	public Map<String, Double> computeFreq(WeightedTagClassJ[] tagClasses) {
		double p = 0.1;
//		println("Initialize p to " + p);

		Set<String> isoformsSet = new HashSet<String>();
		for (int i = 0; i < tagClasses.length; ++i) {
			WeightedTagClassJ tc = tagClasses[i];
			String[] iso = tc.isoforms;
			for (int j = iso.length-1; j >= 0; --j) {
				isoformsSet.add(iso[j]);
			}
		}

		String[] uniqIsoforms = isoformsSet.toArray(new String[isoformsSet.size()]);
		int n = uniqIsoforms.length;

		Map<String, Integer> index = new HashMap<String, Integer>();
		int currentIndex = 0;
		for (int i = 0; i < n; ++i) {
			index.put(uniqIsoforms[i], currentIndex++);
		}

		int[] nrCleaveSites = new int[n];
		double[][] count = new double[n][];
		int maxSites = 0;
		for (int i = 0; i < n; ++i) {
			String isoName = uniqIsoforms[i];
			int nsites = nrCleavageSites.get(isoName);
			if (nsites > maxSites) {
				maxSites = nsites;
			}
			nrCleaveSites[index.get(isoName)] = nsites;
			count[i] = new double[nrCleaveSites[i]];
		}

		double[] f = new double[n];
		Arrays.fill(f, 1.0);
		//		for (i <- uniqIsoforms) f.put(i, 1.0 /*random.uniform(0.0001, 1.0)*/)
		double fSum = Utils.sum(f);

		// speed up future operations by converting isoform names to indices
		for (int i = 0; i < tagClasses.length; ++i) {
			tagClasses[i].indexData(index);
		}

		System.out.print("Iteration: ");
		int steps = 0;
		double maxDiff = 1.0;
		double[] newFreq = new double[n];
		while (maxDiff > 0.00005 && steps < maxSteps) {
			double q = 1 - p;
			double[] qp = preprocessQPowers(q, maxSites);

			// E step
			for (int r = 0; r < tagClasses.length; ++r) {
				WeightedTagClassJ tc = tagClasses[r];

				double sum = 0.0;
				for (int k = tc.size()-1; k >= 0; --k) {
					int i = tc.iso[k];
					int j = tc.pos[k];
					double w = tc.w[k];
					sum +=  f[i] * qp[j] * w; //! make sure j ranges from 0 to nrCleavage(i)-1
				}

				if (sum > EPS) {
					for (int k = tc.size()-1; k >= 0; --k) {
						int i = tc.iso[k];
						int j = tc.pos[k];
						double w = tc.w[k];
						count[i][j] += tc.multiplicity * f[i] * qp[j] * w / sum;
					}
				} else {
					//					println("Give up on " + tc.multiplicity + " tags because sum is " + sum)
				}
			}

			// M step
			double N1 = 0.0;
			double D = 0.0;
			for (int i = 0; i < n; ++i) {
				double r = 1 - qp[nrCleaveSites[i]];
				if (r < EPS) r = 1.0;
				double c = Utils.sum(count[i]); 
				newFreq[i] = c / r;
				D += c;
				N1 += count[i][0];
			}

			// re-estimate p
			p = N1/D;
			System.out.print(steps + " ");
//			println("Number of tags aligning to: first position " + N1
//					+" any position " + D + " p reestimated to " + p);


			// convergence test
			double newFSum = Utils.sum(newFreq);
			maxDiff = 0.0;
			for (int i = 0; i < n; ++i) {
				double oldNormVal = f[i] / fSum;
				if (oldNormVal > EPS) {
					double relDiff = ((newFreq[i] / newFSum) - oldNormVal) / oldNormVal;
					if (relDiff > maxDiff) maxDiff = relDiff;
				}
			}

			// prepare for next iteration
			double[] t = f;
			f = newFreq;
			newFreq = t;

			fSum = newFSum;
			fillWithZeros(count);
			++steps;
			//			println("Step " + steps + " classes " + tagClasses.size)
		}
		if (steps == maxSteps) {
			synchronized(lock) {
				noConverge+=1;
			}
		}
		println("\np converged to " + p);
		Map<String, Double> freq = new HashMap<String, Double>();
		for (int i = 0; i < n; ++i) {
			freq.put(uniqIsoforms[i], f[i]);
		}
		return freq;
	}

	double[] preprocessQPowers(double q, int maxSites) {
		double[] qp = new double[maxSites+1];
		qp[0] = 1;
		for (int i = 1; i <= maxSites; ++i) {
			qp[i] = qp[i-1] * q;
		}
		return qp;
	}

	void fillWithZeros(double[][] m) {
		for (int i = 0; i < m.length; ++i) {
			fillWithZeros(m[i]);
		}
	}

	void fillWithZeros(double[] a) {
		for (int i = 0; i < a.length; ++i) {
			a[i] = 0.0;
		}
	}

}
