#include <iostream>
#include <fstream>

int main(){


  std::fstream outfile;
  outfile.open("lookupList.csv", std::fstream::in | std::fstream::out 
	       | std::fstream::trunc);
  outfile << "ID,events,slopeY,STDY,slopeZ,STDZ\n";
  outfile.close();


}
