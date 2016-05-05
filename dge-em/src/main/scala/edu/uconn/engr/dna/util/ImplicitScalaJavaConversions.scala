package edu.uconn.engr.dna.util

import java.util.HashMap

import edu.uconn.engr.dna.random.RandomPicker

/**
 * User: marius
 * Date: Jul 28, 2010
 */

object ImplicitScalaJavaConversions {
  implicit def myPair2Scala[A, B](p: edu.uconn.engr.dna.util.Pair[A, B]) = scala.Predef.Pair[A, B](p.getFirst, p.getSecond)

  implicit def scala2JavaMap[A, B](m: Map[String, Double]): java.util.Map[java.lang.String, java.lang.Double] = {
    val jm = new HashMap[java.lang.String, java.lang.Double]();
    for ((k, v) <- m) 
    	jm.put(k, v);
    jm
  }

  implicit def randomPickerIntToInteger(p: RandomPicker[Int, Int]) = new RandomPicker[java.lang.Integer, java.lang.Integer] {
    def pick(n: java.lang.Integer): java.lang.Integer = p.pick(n.intValue);
  }
}