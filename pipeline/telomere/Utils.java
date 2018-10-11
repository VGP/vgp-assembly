import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;

import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.BufferedReader;
import java.io.PrintStream;
import java.util.regex.Pattern;
import java.lang.ProcessBuilder;
import java.io.FileInputStream;
import java.io.InputStream;

public class Utils {
   public static final int MBYTES = 1048576;
   public static final int FASTA_LINE_LENGTH = 60;
   public static MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
   public static final String[] FASTA_ENDS = {"trimmed", "contig", "qual", "fasta", "fas", "fna", "fa", "genome"};
   public static final Pattern spaceSplit = Pattern.compile("\\s+");
   public static final Pattern tabSplit = Pattern.compile("\\t+");

   public static class FastaRecord {
      public String header = null;
      public String name = null;
      public String sequence = null;
      public String[] split = null;

      public FastaRecord(String h, String s) {
         split = h.trim().split("\\s+");
         name = split[0];
         header = h;
         sequence = s;
      }
   }

   public static class FastaReader {
      private BufferedReader fileIn = null;
      private String header = "";

      public FastaReader(String file) throws Exception {
         fileIn = Utils.getFile(file, FASTA_ENDS);
      }

      public FastaRecord readRecord() throws Exception {
         String line = null;
         StringBuilder fastaSeq = new StringBuilder();
         FastaRecord result = null;
      
         while ((line = fileIn.readLine()) != null) {
            if (line.startsWith(">")) {
               if (header.length() > 0) {
                  result = new FastaRecord(header, fastaSeq.toString());
                  header = line.substring(1);
                  break;
               }
               header = line.substring(1);
            }
            else {
               fastaSeq.append(line);
            }
         }

         if (result == null && fastaSeq.length() != 0) {
            result = new FastaRecord(header, fastaSeq.toString());
         }

         return result;
      }

      public void close() throws Exception {
        fileIn.close();
      }
   }
   
   public static class Pair {
      public int first;
      public double second;
      public boolean orient;
      public String id;
      
      public Pair(int first, double second) {
         this.first = first;
         this.second = second;
      }
      public Pair(int first, double second, boolean orient) {
         this.first = first;
         this.second = second;
         this.orient = orient;
      }

      public Pair(int first, double second, boolean orient, String id) {
         this.first = first;
         this.second = second;
         this.orient = orient;
         this.id = id;
      }
   }
   
   public enum Translate
   {
       A("T"),
       C("G"),
       G("C"),
       T("A"),
       N("N"),
       R("?"),
       Y("?"),
       K("?"),
       M("?"),
       S("?"),
       W("?"),
       D("?"),
       V("?"),
       H("?"),
       B("?");

       private String other;

       public String getCompliment()
       {
	  if (this == Translate.R) {
	     double rand = Math.random();
             if (rand < 0.5) {
                other = "T";
             }
             else {
                other = "C";
             }
	  }
          else if (this == Translate.Y) {
             double rand = Math.random();
             if (rand < 0.5) {
                other = "G";
             }
             else {
                other = "A";
             }
          }
          else if (this == Translate.K) {
             double rand = Math.random();
             if (rand < 0.5) {
                other = "A";
             }
             else {
                other = "C";
             }
          }
          else if (this == Translate.M) {
             double rand = Math.random();
             if (rand < 0.5) {
                other = "G";
             }
             else {
                other = "T";
             }
          }
          else if (this == Translate.S) {
             double rand = Math.random();
             if (rand < 0.5) {
                other = "C";
             }
             else {
                other = "G";
             }
          }
          else if (this == Translate.W) {
             double rand = Math.random();
             if (rand < 0.5) {
                other = "T";
             }
             else {
                other = "A";
             }
          }
          else if (this == Translate.D) {
             double rand = Math.random()*3;
             if (rand < 1) {
                other = "T";
             }
             else if (rand < 2) {
		other = "C";
             }
             else {
                other = "A";
             }
          }
          else if (this == Translate.V) {
             double rand = Math.random()*3;
             if (rand < 1) {
                other = "T";
             }
             else if (rand < 2) {
                other = "G";
             }
             else {
                other = "C";
             }
          }
          else if (this == Translate.H) {
             double rand = Math.random()*3;
             if (rand < 1) {
                other = "T";
             }
             else if (rand < 2) {
                other = "G";
             }
             else {
                other = "A";
             }
          }
          else if (this == Translate.B) {
             double rand = Math.random()*3;
             if (rand < 1) {
                other = "G";
             }
             else if (rand < 2) {
                other = "C";
             }
             else {
                other = "A";
             }
          }
          return other;
       }
       Translate( String other )
       {
          this.other = other;
       }
   }
   
   public enum ToProtein
   {
      GCT("A"),
      GCC("A"),
      GCA("A"),
      GCG("A"),
      TTA("L"),
      TTG("L"),
      CTT("L"),
      CTC("L"),
      CTA("L"),
      CTG("L"),      
      CGT("R"),
      CGC("R"),
      CGA("R"),
      CGG("R"),
      AGA("R"),
      AGG("R"),
      AAA("K"),
      AAG("K"),
      AAT("N"),
      AAC("N"),
      ATG("M"),
      GAT("D"),
      GAC("D"),
      TTT("F"),
      TTC("F"),
      TGT("C"),
      TGC("C"),
      CCT("P"),
      CCC("P"),
      CCA("P"),
      CCG("P"),
      CAA("Q"),
      CAG("Q"),
      TCT("S"),
      TCC("S"),
      TCA("S"),
      TCG("S"),
      AGT("S"),
      AGC("S"),
      GAA("E"),
      GAG("E"),
      ACT("T"),
      ACC("T"),
      ACA("T"),
      ACG("T"),
      GGT("G"),
      GGC("G"),
      GGA("G"),
      GGG("G"),
      TGG("W"),
      CAT("H"),
      CAC("H"),
      TAT("Y"),
      TAC("Y"),
      ATT("I"),
      ATC("I"),
      ATA("I"),
      GTT("V"),
      GTC("V"),
      GTA("V"),
      GTG("V"),
      TAG("X"),
      TGA("X"),
      TAA("X");
      
      /*
      Ala/A    GCU, GCC, GCA, GCG   
      Leu/L    UUA, UUG, CUU, CUC, CUA, CUG
      Arg/R    CGU, CGC, CGA, CGG, AGA, AGG  
      Lys/K    AAA, AAG
      Asn/N    AAU, AAC    
      Met/M    AUG
      Asp/D    GAU, GAC    
      Phe/F    UUU, UUC
      Cys/C    UGU, UGC    
      Pro/P    CCU, CCC, CCA, CCG
      Gln/Q    CAA, CAG    
      Ser/S    UCU, UCC, UCA, UCG, AGU, AGC
      Glu/E    GAA, GAG    
      Thr/T    ACU, ACC, ACA, ACG
      Gly/G    GGU, GGC, GGA, GGG   
      Trp/W    UGG
      His/H    CAU, CAC    
      Tyr/Y    UAU, UAC
      Ile/I    AUU, AUC, AUA  
      Val/V    GUU, GUC, GUA, GUG
      START    AUG   
      STOP  UAG, UGA, UAA
      */
       private String other;

       public String getProtein()
       {
          return other;
       }
       ToProtein( String other )
       {
          this.other = other;
       }
   }

   public static BufferedReader getFile(String fileName, String postfix) throws Exception {
      String[] arr = new String[1];
      arr[0] = postfix;
      return getFile(fileName, arr);
   }

   public static BufferedReader getFile(String fileName, String[] postfix) throws Exception {
       BufferedReader bf = null;

       bf = new BufferedReader(new InputStreamReader(getFileStream(fileName, postfix, false)));
       bf.ready();
       return bf;
    }

   public static BufferedReader getFile(String fileName, String[] postfix, boolean ignoreUknown) throws Exception {
       BufferedReader bf = null;

       bf = new BufferedReader(new InputStreamReader(getFileStream(fileName, postfix, ignoreUknown)));
       bf.ready();
       return bf;
    }

    public static InputStream getFileStream(String fileName, String[] postfix) throws Exception {
       return getFileStream(fileName, postfix, false);
    }

    public static InputStream getFileStream(String fileName, String[] postfix, boolean ignoreUnknown) throws Exception {
       InputStream is = null;

       if (fileName.endsWith("bz2")) {
          // open file as a pipe
          System.err.println("Running command " + "bzip2 -dc " + new File(fileName).getAbsolutePath() + " |");
          Process p = Runtime.getRuntime().exec("bzip2 -dc " + new File(fileName).getAbsolutePath() + " |");
          is = p.getInputStream();
        } else if (fileName.endsWith("gz")) {
          // open file as a pipe
           System.err.println("Runnning comand " + "gzip -dc " + new File(fileName).getAbsolutePath() + " |");
           Process p = Runtime.getRuntime().exec("gzip -dc " + new File(fileName).getAbsolutePath() + " |");
           is = p.getInputStream();
        } else if (fileName.endsWith("zip")) {
           System.err.println("Runnning comand " + "unzip -p " + new File(fileName).getAbsolutePath() + " |");
           Process p = Runtime.getRuntime().exec("unzip -p " + new File(fileName).getAbsolutePath() + " |");
           is = p.getInputStream();
        } else {
           for (String s : postfix) {
              if (fileName.endsWith(s)){
                 is = new FileInputStream(fileName);
                 return is;
              }
           }
           if (ignoreUnknown == true) {
              is = new FileInputStream(fileName);
              return is;
           }
           System.err.println("Unknown file format " + fileName + " Skipping!");
        }

        return is;
   }

   // add new line breaks every FASTA_LINE_LENGTH characters
   public static String convertToFasta(String supplied) {      
      StringBuffer converted = new StringBuffer();      
      int i = 0;
      String[] split = supplied.trim().split("\\s+");
      if (split.length > 1) { //process as a qual
         int size = 0;
         for (i = 0; i < split.length; i++) {
            converted.append(split[i]);
            size+= split[i].length();
            if (i != (split.length - 1)) {
               if (size >= FASTA_LINE_LENGTH) {
                  size = 0;
                  converted.append("\n");
               } else {
                  converted.append(" ");
               }
            }
         }
      } else { 
         for (i = 0; (i+FASTA_LINE_LENGTH) < supplied.length(); i+= FASTA_LINE_LENGTH) {
            converted.append(supplied.substring(i, i+FASTA_LINE_LENGTH));
            converted.append("\n");         
         }
         converted.append(supplied.substring(i, supplied.length()));
      }      
      return converted.toString();
   }
   
   public static String rc(String supplied) {
      StringBuilder st = new StringBuilder();
      for (int i = supplied.length() - 1; i >= 0; i--) {
         char theChar = supplied.charAt(i);         
         
         if (theChar != '-') {
            Translate t = Translate.valueOf(Character.toString(theChar).toUpperCase());
            st.append(t.getCompliment());
         } else {
            st.append("-");
         }
      }
      return st.toString();
   }

   public static void outputFasta(PrintStream out, String fastaSeq, String ID) {
      outputFasta(out, fastaSeq, null, ID, ">", null, true);
   }

   public static void outputFasta(String fastaSeq, String qualSeq, String ID, String fastaSeparator, String qualSeparator, boolean convert) {
      outputFasta(System.out, fastaSeq, qualSeq, ID, fastaSeparator, qualSeparator, convert);
   }

   public static void outputFasta(PrintStream out, String fastaSeq, String qualSeq, String ID, String fastaSeparator, String qualSeparator, boolean convert) {
      if (fastaSeq.length() == 0) {
         return;
      }

      if (qualSeq != null && qualSeq.length() != fastaSeq.length()) {
         System.err.println("Error length of sequences and fasta for id " + ID + " aren't equal fasta: " + fastaSeq.length() + " qual: " + qualSeq.length());
         System.exit(1);
      }

      out.println(fastaSeparator + ID);
      out.println((convert == true ? Utils.convertToFasta(fastaSeq) : fastaSeq));

      if (qualSeq != null) {
         out.println(qualSeparator + ID);
         out.println((convert == true ? Utils.convertToFasta(qualSeq) : qualSeq));
      }
   }

   public static String getUngappedRead(String fasta) {
      fasta = fasta.replaceAll("N", "");
      fasta = fasta.replaceAll("-", "");
      assert(fasta.length() >= 0);

      return fasta;
   }

   public static int countLetterInRead(String fasta, String letter) {
      return countLetterInRead(fasta, letter, false);
   }

   public static int countLetterInRead(String fasta, String letter, boolean caseSensitive) {
      String ungapped = Utils.getUngappedRead(fasta);
      int len = ungapped.length();
      if (len == 0) { return -1; }
      
      int increment = letter.length();
      int count = 0;

      for (int i = 0; i <= ungapped.length()-increment; i+= increment) {
         if (letter.equals(ungapped.substring(i, i+increment)) && caseSensitive) {
            count++;
         }
         if (letter.equalsIgnoreCase(ungapped.substring(i, i+increment)) && !caseSensitive) {
            count++;
         }
      }
      return count;
   }

   public static String toProtein(String genome, boolean isReversed, int frame) {
      StringBuilder result = new StringBuilder();

      if (isReversed) {
         genome = rc(genome);
      }
      genome = genome.replaceAll("-", "");
      
      for (int i = frame; i < (genome.length() - 3); i += 3) {
         String codon = genome.substring(i, i+3);
         String protein = ToProtein.valueOf(codon).getProtein();
         result.append(protein);
      }
      
      return result.toString();
   }
   
   public static int checkForEnd(String line, int brackets) {
      if (line.startsWith("{")) {
         brackets++;
      }
      if (line.startsWith("}")) {
         brackets--;
      }
      if (brackets == 0) {
         return -1;
      }
      
      return brackets;
   }
   
   public static String getID(String line) {
      String ids[] = line.split(":");
      int commaPos = ids[1].indexOf(",");
      if (commaPos != -1)
    	  return ids[1].substring(1, commaPos).trim();
      else
    	  return ids[1].trim();
   }
   
   public static String getValue(String line, String key) {
      if (line.startsWith(key)) {
         return line.split(":")[1];
      }

      return null;
   }
   
   public static int getOvlSize(int readA, int readB, int ahang, int bhang) {
      if ((ahang <= 0 && bhang >= 0) || (ahang >= 0 && bhang <= 0)) {
         return -1;
      }
      
      if (ahang < 0) {
         return readA - Math.abs(bhang);
      }
      else {
         return readA - ahang;
      }
   }
   
   public static int getRangeOverlap(int startA, int endA, int startB, int endB) {
      int minA = Math.min(startA, endA);
      int minB = Math.min(startB, endB);
      int maxA = Math.max(startA, endA);
      int maxB = Math.max(startB, endB);
      
      int start = Math.max(minA, minB);
      int end = Math.min(maxA, maxB);
      
      return (end-start+1);
   }
}
