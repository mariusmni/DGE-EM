package edu.uconn.engr.dna.dgesim

object TestGeometricProbabilityIndexPicker {
     def main(args: Array[String]) {
    	    val sites = 5
    	    val iterations = 1000000
    	    val count = new Array[Int](5)
    	 	val p = new GeometricProbabilityIndexPicker(0.7, 123, sites);
    	 	for (i <- 0 until iterations) {
    	 		val pos = p.pick(null)
    	 		if (pos >= 0)
    	 			count(pos)+= 1 
    	 	}
    	 	count.foreach(println(_))
    	 	println("Total: " + count.sum)
     }
}