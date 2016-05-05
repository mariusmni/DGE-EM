package edu.uconn.engr.dna.dgeem
import edu.uconn.engr.dna.util.{ ParameterRunnableFactory, BatchThreadPoolExecutor }
import java.io.{ Reader, BufferedReader }
import java.util.{ List, ArrayList }
import java.util.concurrent.ArrayBlockingQueue

class ThreadPooledFastqParser[T](maxThreads: Int,
  queueSize: Int,
  maxBufferSize: Int,
  processorsFactory: ParameterRunnableFactory[List[String], T]) {

  def parse(input: Reader): List[T] = {
    val reader = if (input.isInstanceOf[BufferedReader])
      input.asInstanceOf[BufferedReader]
    else new BufferedReader(input)
    val stringsPerBatch = maxBufferSize - 20
    val queue = new ArrayBlockingQueue[List[String]](queueSize)
    val es = new BatchThreadPoolExecutor[List[String], T](
      maxThreads, queue, processorsFactory, true)
    var l: List[String] = new ArrayList[String](maxBufferSize)
    var tagName = reader.readLine
    while (tagName != null) {
      val tagSequence = reader.readLine
      val plus = reader.readLine
      val qualScores = reader.readLine
      if (l.size >= stringsPerBatch) {
        es.process(l)
        l = es.pollProcessedItemsQueue()
        if (l == null)
          l = new ArrayList[String](maxBufferSize)
        else
          l.clear()
      }
//			l.add(tagName)
      l.add(tagSequence)
	  l.add(qualScores)
      tagName = reader.readLine
    }
    reader.close();
    if (l.size() > 0)
      es.process(l);
    return es.waitForTermination();
  }
}