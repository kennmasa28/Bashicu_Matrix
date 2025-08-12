// C++ headers
#include <cmath>      // sqrt()
#include <csignal>    // ISO C/C++ signal() and sigset_t, sigemptyset() POSIX C extensions
#include <cstdint>    // int64_t
#include <cstdio>     // sscanf()
#include <cstdlib>    // strtol
#include <ctime>      // clock(), CLOCKS_PER_SEC, clock_t
#include <exception>  // exception
#include <ios>        // std::left, std::right
#include <iomanip>    // setprecision()
#include <iostream>   // cout, endl
#include <limits>     // max_digits10
#include <new>        // bad_alloc
#include <string>     // string
#include <sstream>    // ostringstream
#include <fstream>    // ofstream
#include <vector>
#include <algorithm>
#include "athena_arrays.hpp"
#include "BMS.hpp"

int main(int argc, char *argv[]) {
  //--- Step 1. --------------------------------------------------------------------------
  //Read input file
  int n = 1, delta_flag=0, A_flag=0, Bad_flag=0, Good_flag=0, cycle_max = 3, func_flag=2; 
  // Check for command line options and respond.
  for (int i=1; i<argc; i++) {
    // If argv[i] is a 2 character string of the form "-?" then:
    if (*argv[i] == '-'  && *(argv[i]+1) != '\0' && *(argv[i]+2) == '\0') { 
      switch(*(argv[i]+1)) {
        case 'n':
          n = static_cast<int>(std::strtol(argv[++i], nullptr, 10));//argument number          
          break;
        case 'c':
          cycle_max = static_cast<int>(std::strtol(argv[++i], nullptr, 10));//argument number          
          break;
        case 'f':
          func_flag = static_cast<int>(std::strtol(argv[++i], nullptr, 10));//argument number          
          break;
        case 'd':
          delta_flag = 1;        
          break;
        case 'a':
          A_flag = 1;        
          break;
        case 'b':
          Bad_flag = 1;        
          break;
        case 'g':
          Good_flag = 1;        
          break;
        case 'h':
        default:
            std::cout << "Usage: " << argv[0] << " [options] [block/par=value ...]\n";
            std::cout << "Options:" << std::endl;     
            std::cout << "  -n [integer]    substitute argument number into BMS (default: n=1)\n";
            std::cout << "  -c [integer]    decide calculation cycle (default: cycle=3)\n";
            std::cout << "  -f [integer]    Activate function 0:const, 1:n+1, 2:n*n(default)\n";
            std::cout << "  -g              output good part\n";
            std::cout << "  -b              output bad part\n";
            std::cout << "  -d              output Ascension offset delta\n";
            std::cout << "  -a              output Ascension matrix A\n";
            std::cout << "  -h              this help\n";
          return(0);
          break;
      }
    } // else if argv[i] not of form "-?" ignore it here (tested in ModifyFromCmdline)
  }
  

  //--- Step 2. --------------------------------------------------------------------------
  // Set Initial Condition
  Matrix *pm;
  InOut *pio;
  AthenaArray<int> s, g, b, Delta, A, LineAndColumn;
  LineAndColumn.NewAthenaArray(2);
  pio->GetLineAndColumn(LineAndColumn);
  int nline = LineAndColumn(1);
  int ncolumn = LineAndColumn(0);
  s.NewAthenaArray(nline,ncolumn);
  pio->InputMatrix(s);
  std::cout << "\n";
  std::cout << "////// Initial Matrix ////////////////////\n";
  pio->OutputMatrix(s, nline, ncolumn,n);
  
  //--- Step 3. --------------------------------------------------------------------------
  // Calculate BM
  int cycle = 0;
  while (ncolumn > 0 && cycle < cycle_max){
    std::cout << "\n";
    std::cout << "////// step: " << cycle+1 << " ////////////////////\n";
   
    //check the rightest component of s
    int zero_flag = pm->WhetherZero(s, nline, ncolumn);

    if (zero_flag){
      AthenaArray<int> new_s;
      int new_column = ncolumn -1;
      new_s.NewAthenaArray(nline, new_column);
      for (int i=0; i<nline; i++){
        for (int j=0; j<new_column; j++){
          new_s(i,j) = s(i,j);
        }
      }
      ncolumn = new_column;
      n = pm->ActivateFunction(n, func_flag);
      std::cout << "------ New Matrix ------\n";
      pio->OutputMatrix(new_s, nline, ncolumn,n);

      //update
      s.DeleteAthenaArray();
      s.NewAthenaArray(nline, ncolumn);
      for (int i=0; i<nline; i++){
        for (int j=0; j<ncolumn; j++){
          s(i,j) = new_s(i,j);
        }
      }
      new_s.DeleteAthenaArray();
      cycle += 1;

    } else {
    
      //Good Part
      int t = pm->NonZeroLowLine(s,nline,ncolumn);
      int r = pm->DecideBadRoot(s,t,ncolumn);
      if (r==-1){
        std::cout << "Bad root is not found. Finished.\n";
        ncolumn = 0;
        n = pm->ActivateFunction(n, func_flag);
        break;
      }
      g.NewAthenaArray(nline,r);
      pm->DecideGoodPart(g, s, r, nline);
      
      //Delta
      Delta.NewAthenaArray(nline,1);
      pm->DecideDelta(Delta, s, nline, ncolumn, t, r);

      //A
      A.NewAthenaArray(nline, ncolumn-r-1);
      pm->DecideA(A, s, nline, ncolumn, t, r);
      
      //Bad Part
      int fn = pm->ActivateFunction(n, func_flag);
      b.NewAthenaArray(fn+1 ,nline, ncolumn-r-1);
      pm->DecideBadPart(b,s, Delta, A, nline, ncolumn-r-1, r, fn);
      
      //New BM and new n
      AthenaArray<int> new_s;
      int new_column = r + (ncolumn - r -1 )*(fn+1);
      new_s.NewAthenaArray(nline, new_column);
      pm->NewMatrix(new_s, g, b, nline, ncolumn-r-1, r, fn);
      n = fn;

      //Output
      if (Good_flag){
        std::cout << "------ Good Part ------\n";
        pio->OutputMatrix(g, nline, r,-1);
      }
      if (Bad_flag)
      {
        std::cout << "------ Bad Part ------\n";
        pio->Outputb(b, 0, nline, ncolumn-r-1);
      }
      if (delta_flag){
        std::cout << "------ Ascension offset delta ------\n";
        pio->OutputMatrix(Delta, nline, 1,-1);
      }
      if (A_flag){
        std::cout << "------ Ascension matrix A ------\n";
        pio->OutputMatrix(A, nline, ncolumn-r-1,-1);
      }
      std::cout << "------ New Matrix ------\n";
      pio->OutputMatrix(new_s, nline, new_column,n);
      
      //update
      ncolumn = new_column;
      s.DeleteAthenaArray();
      b.DeleteAthenaArray();
      g.DeleteAthenaArray();
      A.DeleteAthenaArray();
      Delta.DeleteAthenaArray();
      s.NewAthenaArray(nline, ncolumn);
      for (int i=0; i<nline; i++){
        for (int j=0; j<ncolumn; j++){
          s(i,j) = new_s(i,j);
        }
      }
      new_s.DeleteAthenaArray();
      n = fn;
      cycle += 1;
    }//if...else
  }//while

  if (ncolumn==0){
    std::cout << "\n";
    std::cout << "////// Calculation is finished !! ////////////////////\n";
    std::cout << "final output = " << n << "\n";

  }else if (cycle == cycle_max){
    std::cout << "\n";
    std::cout << "////// Calculation cycle is up to limit ////////////////////\n";
  }
  
}

std::vector<std::string> InOut::split(std::string str, std::string separator) {
  if (separator == "") return {str};
  std::vector<std::string> result;
  std::string tstr = str + separator;
  unsigned long l = tstr.length(), sl = separator.length();
  std::string::size_type pos = 0, prev = 0;
  
  for (;pos < l && (pos = tstr.find(separator, pos)) != std::string::npos; prev = (pos += sl)) {
    result.emplace_back(tstr, prev, pos - prev);
  }
  
  return result;
}

void InOut::GetLineAndColumn( AthenaArray<int> &LineAndColumn){
  std::ifstream ifs("BMS.txt");
  
  if (ifs.fail()) {
    std::cerr << "Failed to open the input file" << std::endl;
    return;
  }
  std::string str;
  int j = 0;

  while (std::getline(ifs, str)) {
    //全文ifsの一行ずつを収めたstrを,の前後で分割
    std::vector<std::string> ary = split(str,",");
    // remove spaces
    LineAndColumn(0) = ary.size();
    for (unsigned int i = 0; i < ary.size(); i++) {
        std::string temp = ary[i];
        temp.erase(remove(temp.begin(), temp.end(),' '), temp.end());
    }
    j+=1;

  }
  LineAndColumn(1) = j; 
  return;
}

void InOut::InputMatrix(AthenaArray<int> &s){
  
  std::ifstream ifs("BMS.txt");
  
  if (ifs.fail()) {
      std::cerr << "Failed to open the input file" << std::endl;
      return;
  }
  std::string str;
  int j = 0;

  while (std::getline(ifs, str)) {
    //全文ifsの一行ずつを収めたstrを,の前後で分割
    std::vector<std::string> ary = split(str,",");
    // remove spaces
    for (unsigned int i = 0; i < ary.size(); i++) {
        std::string temp = ary[i];
        temp.erase(remove(temp.begin(), temp.end(),' '), temp.end());
        s(j,i) = stoi(temp);
    }
    j+=1;
  }
  return;
}

void InOut::OutputMatrix(AthenaArray<int> &s, int nline, int ncolumn, int n){
  for (int i=0; i<nline; i++){
    for (int j=0; j<ncolumn; j++){
      std::cout << std::right << std::setw(3) << s(i,j);
      if (j<ncolumn-1) std::cout  << ",";
    }
    std::cout << "\n";
  }
  if (n>0) std::cout << "[" << n << "]\n";
}

void InOut::Outputb(AthenaArray<int> &b, int a, int nline, int maxcolumn){
  for (int i=0; i<nline; i++){
    for (int j=0; j<maxcolumn; j++){
      std::cout << std::right << std::setw(3) << b(a,i,j);
      if (j<maxcolumn-1) std::cout << ",";
    }
    std::cout << "\n";
  }
}

int Matrix::WhetherZero(AthenaArray<int> &s, int nline, int ncolumn){
  for (int i=0; i<nline; i++){
    if (s(i,ncolumn-1)!=0) {
      return 0;
    }
  }
  return 1;
}

int Matrix::NonZeroLowLine(AthenaArray<int> &s, int nline, int ncolumn){
  for (int i=0; i<nline; i++){
    int t = nline - i - 1;
    if (s(t,ncolumn-1)!=0) {
      return t;
    }
  }
  return -1;
}

int Matrix::DecideBadRoot(AthenaArray<int> &s, int t, int ncolumn){
  return ColumnOfParent(s, t, ncolumn-1);
}

int Matrix::ColumnOfParent(AthenaArray<int> &s, int x, int y){
  if (x==0 && y>0){
    for (int i=0; i<y; i++){
      int r = y - i - 1;
      if (s(x,r) < s(x,y)){
        //std::cout << "r = " << r << "\n";
        return r;
      }
    }
    return -1;

  } else if (x>0 && y>0){
    for (int i=0; i<y; i++){
      int r = y - i - 1;
      if (s(x,r) < s(x,y) && IsColumnOfDirectAncestor(s, r, x-1, y)){
        //std::cout << "r = " << r << "\n";
        return r;
      }
    }
    return -1;
  }
}

int Matrix::IsColumnOfDirectAncestor(AthenaArray<int> &s, int column, int x, int y){
  //Judge "is column r direct ancestor of s(line=x, column=y)?"

  int j=y;
  while(j>=0){
    if(column==j) return 1;
    if(column==ColumnOfParent(s, x, j)) return 1;
    j = ColumnOfParent(s, x, j);
  }
    
  return 0;
}

void Matrix::DecideGoodPart(AthenaArray<int> &g, AthenaArray<int> &s, int r, int nline){
  for (int i=0; i<nline; i++){
    for (int j=0; j<r; j++){
      g(i,j) = s(i,j);
    }
  }
}

void Matrix::DecideDelta(AthenaArray<int> &Delta, AthenaArray<int> &s, int nline, int ncolumn, int t, int r){
  for (int i=0; i<nline; i++){
    if (i<t) Delta(i) = s(i, ncolumn-1) - s(i, r);
    if (i>=t) Delta(i) = 0;
    //std::cout << Delta(i) << "\n";
  }
}

void Matrix::DecideA(AthenaArray<int> &A, AthenaArray<int> &s, int nline, int ncolumn, int t, int r){
  for (int i=0; i<nline; i++){
    for (int j=0; j<ncolumn - r - 1; j++){
      // std::cout << "r=" << r << "| IsCD=" << IsColumnOfDirectAncestor(s, r, i, j+r)<< "\n";
      if (IsColumnOfDirectAncestor(s, r, i, j+r)) A(i,j) = 1;
      if (!IsColumnOfDirectAncestor(s, r, i, j+r)) A(i,j) = 0;
    }
  }
}

void Matrix::DecideBadPart(AthenaArray<int> &b , AthenaArray<int> &s,  AthenaArray<int> &Delta, AthenaArray<int> &A, int nline, int maxcolumn, int r, int fn){
  for (int a=0; a<=fn; a++){
    for (int i=0; i<nline; i++){
      for (int j=0; j<maxcolumn; j++){
        b(a,i,j) = s(i,r+j) + a*Delta(i)*A(i,j);
      }
    }
  }
}

int Matrix::ActivateFunction(int n, int func_flag){
  switch (func_flag)
  {
    case 0:
      return n;
      break;
    case 1:
      return n+1;
      break;
    case 2:
      return n*n;
      break;
    default:
      break;
  }
  return n*n;
}

void Matrix::NewMatrix(AthenaArray<int> &new_s,  AthenaArray<int> &g, AthenaArray<int> &b, int nline, int maxcolumn, int r, int fn){
  //std:: cout << nline << " " << r << "\n";
  for (int i0=0; i0<nline; i0++){
    for (int j0=0; j0<r; j0++){
      new_s(i0,j0) = g(i0,j0);
    }
  }

  for (int i=0; i<nline; i++){
    for (int a=0; a<=fn; a++){
      for (int j=0; j<maxcolumn; j++){
        new_s(i,r+j+maxcolumn*a) = b(a,i,j);
        //if (i==1) std::cout << new_s(i,r+j+maxcolumn*a) << "\n";
      }
    }
  }

}