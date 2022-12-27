#include <iostream>       // std::cout
#include <string>         // std::string
#include <fstream>

std::string rc(std::string str) {
  std::string DNAseq(str);

  size_t c = 0;
  // reverse
  for(int i = str.length()-1; i>=0; i--) {
     DNAseq[c++] = str[i];
  }
  // complement
  for (std::size_t i = 0; i < DNAseq.length(); ++i){
        switch (DNAseq[i]){
        case 'A':
            DNAseq[i] = 'T';
            break;
        case 'C':
            DNAseq[i] = 'G';
            break;
        case 'G':
            DNAseq[i] = 'C';
            break;
        case 'T':
            DNAseq[i] = 'A';
            break;
        }
    }
    return DNAseq;
}

void find (std::string str, std::string name, std::string str2) {
  std::size_t pos = str.find(str2);
  while (pos != std::string::npos) {
     std::size_t len = 0;
     std::cout << name << "\t" << str.length() << "\t0\t" << pos;
     while (str.substr(pos, str2.length()) == str2) {
        len+=str2.length();
        pos+=str2.length();
     }
     std::cout << "\t" << pos << "\t" << len<< std::endl;
     pos = str.find(str2, pos+1);
  }

  // now ref search
  std::string rev = rc(str2);
  pos = 0;
  pos = str.find(rev);
  while (pos != std::string::npos) {
     std::size_t len = 0;
     std::cout << name << "\t" << str.length() << "\t1\t" << pos;
     while (str.substr(pos, rev.length()) == rev) {
        len+=rev.length();
        pos+=rev.length();
     }
     std::cout << "\t" << pos << "\t" << len << std::endl;
     pos = str.find(rev, pos+1);
  }
}

int main(int argc, char * argv[])
{
  std::string str ("");
  std::ifstream infile(argv[1]);
  std::string line;
  std::string header;

  if (argc < 2) {
     std::cerr << "Error: invalid number of parameters" << std::endl;
     std::cerr << "Usage: find <input fasta> [optional sequence to search for, default is vertebrate TTAGG" << std::endl;
     exit(1);
   }

  while (std::getline(infile, line)) {
     if (line.find(">") == std::string::npos) {
        str.append(line);
     } else {
       if (str.length() > 0) { find(str, header, (argc >=3 ? argv[2] : "TTAGGG")); }
       str = "";
       header=line;
     }
  }
  if (str.length() > 0) { find(str, header, (argc >= 3 ? argv[2] : "TTAGGG")); }
  infile.close();

  return 0;
}
