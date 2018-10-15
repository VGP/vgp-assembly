#!/bin/bash

javac *.java
jar cf telomere.jar *.class
rm *.class

g++ find_telomere.c -o find_telomere
