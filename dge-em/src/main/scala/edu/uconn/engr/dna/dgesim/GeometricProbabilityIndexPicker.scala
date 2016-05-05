package edu.uconn.engr.dna.dgesim

import edu.uconn.engr.dna.random.RandomPicker
import util.Random
/**
 * Last site has probability p,
 * Previous has p*(1-p)
 * Previous has p*(1-p)^2
 * ...
 * User: marius
 * Date: Jul 28, 2010
 */

class GeometricProbabilityIndexPicker(p: Double, seed: Int, numberOfSites:Int) extends RandomPicker[Int, Null] {
  import scala.math._
  import GeometricProbabilityIndexPicker._
  val logOneMinusP = log(1-p)
  val random = new Random(seed)

  def pick(dummy: Null): Int = {
    val r = random.nextDouble

    if (r < EPS) // very small value
      if (r <= power(1 - p, numberOfSites))
        return -1

    val k = floor(log(r) / logOneMinusP).intValue
    if (k >= numberOfSites)
      -1
    else
      k
  }
}

object GeometricProbabilityIndexPicker {
  val EPS = 0.000000001
  
  /**
   * Log time exponentiation
   */
  def power(a: Double, n: Int): Double = {
    if (n == 0)
      1
    else {
      val b = power(a, n / 2)
      if (n % 2 == 1)
        b * b * a
      else
        b * b
    }
  }
}
