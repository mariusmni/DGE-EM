package edu.uconn.engr.dna.dgesim

import java.io.InputStream
import java.util.Properties
import edu.uconn.engr.dna.probability.{GeometricProbabilityDistribution, UniformProbabilityDistribution, ProbabilityDistribution}

/**
 * User: marius
 * Date: Jul 26, 2010
 */

class SimulationParams(inputStream: InputStream) extends Properties {
  load(inputStream)

  override def getProperty(s: String): String = {
	  val v = super.getProperty(s)
	  if (v == null) null else v.trim
  }
  
  def getGTF: String = getProperty("gtf")

  def getClusters: String = getProperty("clusters")

  def isoformDistribution: ProbabilityDistribution = {
    val parts = getProperty("isoformDistribution").split(",")
    val name = parts(0).trim
    name match {
      case "uniform" => new UniformProbabilityDistribution(randomNumberGeneratorSeed)
      case "geometric" => {
        var p = 0.5
        if (parts.length > 0)
          p = parts(1).trim.toDouble
        new GeometricProbabilityDistribution(p, randomNumberGeneratorSeed)
      }
    }
  }

  def randomNumberGeneratorSeed: Int = getProperty("randomNumberGeneratorSeed").toInt

  def clusterDistributionFileName: String = getProperty("clusterDistributionFileName")

  def isoformSequencesFile: String = getProperty("isoformSequencesFile")

  def tagLength: Int = getProperty("tagLength").toInt

  def numberOfTags: Int = getProperty("numberOfTags").toInt

  def p: Double = getProperty("p").toDouble

  def outputIsoformWeights: Boolean = getProperty("outputIsoformWeights", "true").toBoolean

  def cleavePatterns: List[String] = getProperty("cleavePatterns").trim.split(",").toList

  def printToConsole: Boolean = {
	  val p = getProperty("printToConsole")
	  if (p == null) false else p.toBoolean
  }

}