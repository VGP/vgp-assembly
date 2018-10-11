#!/bin/bash

javac *.java
jar cf telomere.jar *.class
rm *.class

g++ find.c -o find
