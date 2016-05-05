@echo off
java -cp "$INSTALL_PATH/lib/$jarname" edu.uconn.engr.dna.format.GenerateKnownGeneSequences %*
