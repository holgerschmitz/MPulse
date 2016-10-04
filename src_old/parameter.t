
//the rebuild method template. called by Rebuildable objects
//gets *pValue from stream returns next token
template<class Type>
std::string ParameterValue<Type>::Rebuild (std::istream& in) {
  in >> *pValue;
  std::string strToken;
  in >> strToken;
  return strToken;
}


template<class Type, class BaseType>
std::string ParameterRebuild<Type,BaseType>::Rebuild (std::istream& in) {
  BaseType *val = NewInstance();
  if (value) (*value) = val;
  else if (values) values->push_back(val);
  

  // Task aus der Datei wiederherstellen
  std::string strToken;
  in >> strToken;
  if ((strToken != "{") || in.eof()) return strToken;
  std::string nextToken = val->Rebuild(in);
  val->finalize();
  return nextToken;
}
