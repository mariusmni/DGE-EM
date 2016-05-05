package edu.uconn.engr.dna.annotation;

public @interface SimpleAnnotationExample {
	int    id();
	String synopsis();
	String engineer() default "[unassigned]"; 
	String date() default "[unimplemented]"; 
}
