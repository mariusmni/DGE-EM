package edu.uconn.engr.dna.format;

import edu.uconn.engr.dna.util.AbstractRandomAccessMap;



public class Isoforms extends AbstractRandomAccessMap<String, Isoform, Isoform> {

	@Override
	protected Isoform addToExistingGroupOrCreateNew(String key, Isoform isoform, Isoform value) {
		return value;
	}

	public Iterable<Isoform> isoformIterator() {
		return groupIterator();
	}

}
