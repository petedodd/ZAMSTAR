//#######################
#include "person.h"
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//HOUSEHOLD MEMBER FUNS
int household::add( unsigned int who, unsigned int sex ){
  //  if ( parz.hhFLAG ){
    members.push_back( who );
    //  }
  int j = 0;
  switch( sex ){
  case 0://woman
    nwomen++;
    j++;
    break;
  case 1:
    nmen++;
    j++;
    break;
  default:
    j = -1;
    break;
  }
  return j;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int household::remove( unsigned int who, unsigned int sex){
  //  if ( parz.hhFLAG ){//don't want to spend ages looking for them when hhs don't exist...
    members.remove( who );
    //  }
  int j = 0;
  switch( sex ){
  case 0://woman
    nwomen--;
    j++;
    break;
  case 1:
    nmen--;
    j++;
    break;
  default:
    j = -1;
    break;
  }
  return j;
}
//end household member funs
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

