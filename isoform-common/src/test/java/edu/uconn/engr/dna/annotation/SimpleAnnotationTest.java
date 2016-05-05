package edu.uconn.engr.dna.annotation;

import static junit.framework.Assert.assertEquals;

import java.util.ArrayList;

import org.junit.Test;

public class SimpleAnnotationTest {

	@SimpleAnnotationExample (
			id = 1234,
			engineer = "marius",
			date = "6/6/10",
			synopsis = "annotation test"
	)
	public void olala() {

	}

	@Test(expected=IndexOutOfBoundsException.class, timeout=100)
	public void testSum() {
		new ArrayList<String>().get(0);
		assertEquals(2, 2);
	}
}
