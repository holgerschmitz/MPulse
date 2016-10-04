#include "rebuild.h"
#include <algorithm>
//for the count algorithm
//-----------------------------------------------------------------------------
ParameterMap* Rebuildable::MakeParamMap (ParameterMap* pm) {
  //if not present, allocate
  if (NULL == pm) pm = new ParameterMap;
  return pm;
}
//ParameterMap is a std::map<std::string,WParameter> as found in parameter.h, while WParameter 
// is a wrapped pointer to parameter
// for ParameterMap see parameter.h
//-----------------------------------------------------------------------------
//Rebuild
std::string Rebuildable::Rebuild (std::istream& in) {
  //set up parameter map
  ParameterMap* pm = MakeParamMap();
  // get first token 
  std::string strToken;
  in >> strToken;
  //outer while loop
  //get next token until EOF or "}", remove commentary lines
  while ((strToken != "}") && (!in.eof())) {
    //-------------------------------------------------------------------------
    // loop for removing commentaries from stream
    while (strToken == "//") { 
      char ch;
      do {
		    in.get(ch);
      } while ((ch != '\n') && (!in.eof()));
      //get next token unless EOF is encountered
      if (!in.eof()) in >> strToken;
    }
    //end of commentary loop
    //-----------------------------------------------------------------------
    //test on not (end of block)
    if ("}" != strToken) {
      //----------------------------------------------------------------------
      //handling errors and unexpected tokens
      //if its not found in the parameter map, write it to strParam for displaying error messages
      if (0 == pm->count(strToken)) { 
        std::string strParam = strToken;
        		//get next token
				in >> strToken;
				//if a "{" token is encountered, skip the contents
				// since this is a task, it should be skipped (?)
				if (strToken == "{") { 
				  std::cerr << "Unknown Task " << strParam << ". Skipping... " << std::endl;
				  int nLevel = 0;
				  do {
				    in >> strToken;
				    if (strToken == "{") nLevel++;
				    if (strToken == "}") nLevel--;
				  } while (((strToken != "}") || (nLevel >= 0)) && (!in.eof()));
		  		  //if EOF encountered, display error message, else get next token
          		  if (in.eof()) 
				    std::cerr << "Unexpected end of file." << std::endl;
				  else
				    in >> strToken;
				}
				else 
				  // Display name of parameter if it is not found in the map
				  // But do not break, since a correct one could follow
				  std::cerr << "Unknown Parameter " << strParam << std::endl;
      }
      //-----------------------------------------------------------------------
      // here the real rebuilding is done
      // this else belongs to if (0 == pm->count(strToken))!
      // if all is as expected, write data to the parameters 
	// by calling its rebuild method
      else {
        // get the parameter object from the map via index operator
        // dereference pm[token] --> get wrapped pointer to parameter 
	  //--> cast to normal pointer to parameter 
	  // --> assign to new pointer to parameter 
        // c-style cast! (reinterpret_cast<*parameter>(*pm) ?)
        Parameter *par = (Parameter*)(*pm)[strToken];
        // deference and call the rebuild method  of this object with the "in" filestream
	  // parameter's rebuild member gets the next token from the file stream, 
	  // until a new block is encountered, see also parameter.h
        strToken = par->Rebuild(in);
      }
      //-----------------------------------------------------------------------
    }
  }
  //---------------------------------------------------------------------------
  // end of outer while loop
  // if end of block reached, get next token
  if (strToken == "}") in >> strToken;
  //remove commentary from instream
  while (strToken == "//") { 
    char ch;
    do {
      in.get(ch);
    } while ((ch != '\n') && (!in.eof()));
    //get next token until EOF
    if (!in.eof()) 
      in >> strToken;
    else
      break;
  }
  //---------------------------------------------------------------------------
  // clean up
  // clean pm
  while (!pm->empty()) {
    ParameterMap::iterator iter = pm->begin();
    //c-style cast to normal, non-wrapped pointer to parameter
    delete ((Parameter*)(*iter).second);
    pm->erase(iter);
  }
  //deallocate pm
  delete pm;
  //---------------------------------------------------------------------------
  //return token
  return strToken;
}



