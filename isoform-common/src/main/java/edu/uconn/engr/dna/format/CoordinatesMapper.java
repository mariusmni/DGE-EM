package edu.uconn.engr.dna.format;

import edu.uconn.engr.dna.util.IntervalVisitor;
import edu.uconn.engr.dna.util.Intervals;

public class CoordinatesMapper {


	/**
	 * <p>
	 * Converts a genome position to a position in the coordinates of a isoform.
	 * The genome position must actually fall within the positions covered by the
	 * isoform, otherwise -1 is returned. 
	 * </p>
	 * 
	 * <p>
	 * For example if the isoform covers intervals [2..5] and [8..12] in 
	 * the genome, then coordinate 2 would be the 1st position in the isoform, 
	 * coordinate 4 would be the 3rd, coordinate 6 would not fall in the isoform
	 * and coordinate 8 would be the 5th position covered by the isoform. In other 
	 * words, for genome positions 2, 4, 6, 8 and the isoform given above, this method
	 * would return 1, 3, -1, and 5 respectively.
	 * </p>
	 * 
	 * @param genomePosition a coordinate in the genome
	 * @param isoformIntervalsInGenome the intervals covered by the isoform in
	 * the genome
	 * @return the number of the genome position in isoform coordinates (1 based), or -1
	 * if the genome position is not covered by the isoform
	 */
	public static int genomeToIsoformCoordinates(int genomePosition, Intervals isoformIntervalsInGenome) {
		int pos = 1;
		int n = isoformIntervalsInGenome.size();
		for (int i = 0; i < n; ++i) {
			int intStart = isoformIntervalsInGenome.getStart(i);
			if (intStart > genomePosition) {
				return -1;
			}
			int intEnd = isoformIntervalsInGenome.getEnd(i);
			if (intEnd >= genomePosition) {
				return pos +  genomePosition - intStart;
			}
			pos += intEnd - intStart + 1;
		}
		return -1;
	}

	/**
	 * Converts an interval from isoform coordinates to genome coordinates
	 * 
	 * @param isoformIntervalsInGenome the intervals covered by the isoform in the genome
	 * @param readStart the start of the subinterval in the isoform that we are interested 
	 * to convert to genome coordinates
	 * @param readEnd the end of the subinterval in the isoform that we are interested 
	 * to convert to genome coordinates
	 * @return a set of Intervals representing the genome intervals covered by the isoform
	 * subinterval
	 */
	public static Intervals isoformToGenomeCoordinates(Intervals isoformIntervalsInGenome, int readStart, int readEnd) {
		final Intervals intervals = new Intervals();
		visitGenomeIntervals(isoformIntervalsInGenome, readStart, readEnd, new IntervalVisitor() {
			
			@Override
			public void visit(int isoformIntervalStart, int isoformIntervalEnd,
					int genomeIntervalStart, int genomeIntervalEnd) {
//				System.out.printf("{%d %d} ", isoformIntervalStart, isoformIntervalEnd);
				intervals.add(genomeIntervalStart, genomeIntervalEnd);
			}
		});
//		System.out.println();
		return intervals;
	}

	public static Intervals isoformToGenomeCoordinatesOld(Intervals isoformIntervalsInGenome, int readStart, int readEnd) {
		if (readStart <= 0) {
			throw new IllegalArgumentException("Read start must be greater than 0");
		}
		int sum = 0;
		int p = 0;
		while (sum < readStart) {
			sum += isoformIntervalsInGenome.length(p++);
		}
		int sump = sum - isoformIntervalsInGenome.length(p-1);
		
		int q = p;
		while (sum < readEnd) {
			sum += isoformIntervalsInGenome.length(q++);
		}
		
		int sumq = sum - isoformIntervalsInGenome.length(q-1);
		--p;
		--q;
		
		if (p == q) {
			return new Intervals(new int[]{isoformIntervalsInGenome.getStart(p) + readStart - sump - 1},
								 new int[]{isoformIntervalsInGenome.getStart(p) + readEnd - sump - 1});
		} else {
			Intervals output = new Intervals();
			output.add(isoformIntervalsInGenome.getStart(p) + readStart - sump - 1, isoformIntervalsInGenome.getEnd(p));
			for (int i = p+1; i < q; ++i) {
				output.add(isoformIntervalsInGenome.getStart(i), isoformIntervalsInGenome.getEnd(i));
			}
			output.add(isoformIntervalsInGenome.getStart(q), isoformIntervalsInGenome.getStart(q) + readEnd - sumq - 1);
			return output;
		}
	}

	public static int isoformToGenomeCoordinate(Intervals isoformIntervalsInGenome, int readStart) {
		if (readStart <= 0) {
			throw new IllegalArgumentException("Position start must be greater than 0");
		}
		int sum = 0;
		int p = 0;
		while ((sum += isoformIntervalsInGenome.length(p++)) < readStart) {
		}
		int sump = sum - isoformIntervalsInGenome.length(--p);
		return isoformIntervalsInGenome.getStart(p) + readStart - sump - 1;
	}

	/**
	 * Converts the isoform subinterval [readStart, readEnd] into a set of genome
	 * intervals and calls the visit method on the IntervalVisitor for each of the
	 * resulting intervals.
	 *    
	 * @param isoformIntervalsInGenome the intervals covered by the isoform in the genome
	 * @param readStart the start of the subinterval in the isoform that we are interested 
	 * to convert to genome coordinates
	 * @param readEnd the end of the subinterval in the isoform that we are interested 
	 * to convert to genome coordinates
	 * @param iv an IntervalVisitor which will be notified of each genome interval
	 * resulted from the conversion
	 */
	public static void visitGenomeIntervals(Intervals isoformIntervalsInGenome, 
			int readStart, int readEnd, IntervalVisitor iv) {
		if (readStart <= 0) {
			throw new IllegalArgumentException("Read start must be greater than 0");
		}
		int sum = 0;
		int p = 0;
		while ((sum += isoformIntervalsInGenome.length(p++)) < readStart) {
		}
		int sump = sum - isoformIntervalsInGenome.length(p-1);
		
		int q = p;
		while (sum < readEnd) {
			sum += isoformIntervalsInGenome.length(q++);
		}
		
		int sumq = sum - isoformIntervalsInGenome.length(q-1);
		--p;
		--q;
		
		if (p == q) {
			int genostart = isoformIntervalsInGenome.getStart(p) + readStart - sump - 1;
			int genomeend = isoformIntervalsInGenome.getStart(p) + readEnd - sump - 1;
			iv.visit(readStart, readEnd, 
					genostart, genomeend);
		} else {
			int readPos = readStart;
			readPos = visit(readPos, 
					isoformIntervalsInGenome.getStart(p) + readStart - sump - 1, 
					isoformIntervalsInGenome.getEnd(p), 
					iv);
			for (int i = p+1; i < q; ++i) {
				readPos = visit(readPos, 
						isoformIntervalsInGenome.getStart(i), 
						isoformIntervalsInGenome.getEnd(i), 
						iv);
			}
			readPos = visit(readPos, 
					isoformIntervalsInGenome.getStart(q), 
					isoformIntervalsInGenome.getStart(q) + readEnd - sumq - 1, 
					iv);
		}

	}

	private static int visit(int readPos, int genostart, int genomeend, IntervalVisitor iv) {
		int readPosEnd = readPos + genomeend - genostart;
		iv.visit(readPos, readPosEnd, genostart, genomeend);
		return readPosEnd + 1;
	}

	/**
	 * Converts a set of intervals from isoform coordinates to genome doordinates
	 * 
	 * @param isoformIntervalsInGenome the intervals covered by the isoform in the genome
	 * @param isoformSubIntervals the intervals in the isoform that we are interested to convert
	 * to genome coordinates
	 * @return a set of intervals obtained by translating every interval from isoformSubIntervals
	 * into genome coordinates, and reuniting all the resulted intervals
	 */
	public static Intervals isoformToGenomeCoordinates(Intervals isoformIntervalsInGenome, Intervals isoformSubIntervals) {
		Intervals result = new Intervals();
		for (int i = 0; i < isoformSubIntervals.size(); ++i) {
			Intervals genomeIntervalsForSegment = isoformToGenomeCoordinatesOld(isoformIntervalsInGenome, 
					isoformSubIntervals.getStart(i), isoformSubIntervals.getEnd(i));
			result.append(genomeIntervalsForSegment);
		}
		return result;
	}

	public static String createGenomeCigar(Intervals intervals) {
		StringBuilder cigar = new StringBuilder();
		cigar.append(intervals.length(0));
		cigar.append('M');
		for (int j = 1; j < intervals.size(); ++j) {
			cigar.append(intervals.getGapSize(j));
			cigar.append('N');
			cigar.append(intervals.length(j));
			cigar.append('M');
		}

		return cigar.toString();
	}
	
}
