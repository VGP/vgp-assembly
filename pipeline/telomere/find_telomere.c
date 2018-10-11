#include <iostream>       // std::cout
#include <string>         // std::string
#include <fstream>

void find (std::string str, std::string name) {
  std::string str2 ("TTAGGG");

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
  str2 = "CCCTAA";
  pos = 0;
  pos = str.find(str2);
  while (pos != std::string::npos) {
     std::size_t len = 0;
     std::cout << name << "\t" << str.length() << "\t1\t" << pos;
     while (str.substr(pos, str2.length()) == str2) {
        len+=str2.length();
        pos+=str2.length();
     }
     std::cout << "\t" << pos << "\t" <<  len << std::endl;
     pos = str.find(str2, pos+1);
  }
}

int main(int argc, char * argv[])
{
  std::string str ("");
  std::ifstream infile(argv[1]);
  std::string line;
  std::string header;
  while (std::getline(infile, line)) {
     if (line.find(">") == std::string::npos) {
        str.append(line);
     } else {
       if (str.length() > 0) { find(str, header); }
       str = "";
       header=line;
     }
  }
  if (str.length() > 0) { find(str, header); }
  infile.close();

  return 0;
}
