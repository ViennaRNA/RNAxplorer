/**********************************************/
/* BEGIN interface for distorted_sampling functions   */
/**********************************************/

//no default constructor
//%nodefaultdtor gridLandscapeT;

/*
%extend gridLandscapeT {
	~gridLandscapeT(){
	  	for(int i = 0; i < $self->size1; i++){
	  		free($self->landscape[i]);
	  	}
		free($self->landscape);		
	  	free($self);
	}
}
*/

%typemap(out) gridLandscapeT* estimate_landscape %{
	/*
	int ownership = 1; // 1 == python.
	gridLandscapeT* gp = $1;
	PyObject* res = SWIG_NewPointerObj(gp, SWIGTYPE_p_gridLandscapeT, ownership);
	$result = res;
	*/
	
	//get number of cells
	int numberOfCells = 0;
	for(int i = 0; i < $1->size1;i++){
		for(int j = 0; j < $1->size2;j++){
			if($1->landscape[i][j].num_structs > 0){
				numberOfCells++;
			}
		}
	}
	
	int cellIndex = 0;
	PyObject* res = PyList_New(numberOfCells);
	for(int i = 0; i < $1->size1; i++){
		for(int j = 0; j < $1->size2;j++){
			if($1->landscape[i][j].num_structs > 0){
				PyObject* cellList = PyTuple_New(4);
				
				PyObject* pK = PyInt_FromLong($1->landscape[i][j].k);
				PyObject* pL = PyInt_FromLong($1->landscape[i][j].l);
				
				double num = $1->landscape[i][j].mfe;
				PyObject* pMFE = PyFloat_FromDouble(num);
				
				PyObject* structureList = PyList_New($1->landscape[i][j].num_structs);
				for(int k = 0; k < $1->landscape[i][j].num_structs; k++){
					PyObject* pStr = PyString_FromString($1->landscape[i][j].structures[k]);
					PyList_SetItem(structureList,k,pStr);
					free($1->landscape[i][j].structures[k]);
				}
				
				PyTuple_SET_ITEM(cellList,0,pK);
				PyTuple_SET_ITEM(cellList,1,pL);
				PyTuple_SET_ITEM(cellList,2,pMFE);
				PyTuple_SET_ITEM(cellList,3,structureList);
				PyList_SetItem(res,cellIndex,cellList);
				cellIndex++;
			}
			free($1->landscape[i][j].structures);
		}
		free($1->landscape[i]);
	}
	free($1->landscape);
	free($1);
	$result=res;
%}

/*
%typemap(freearg) char ** {
  free($1);
}

%typemap(out) char ** {
  int len,i;
  len = 0;
  while ($1[len]) len++;
  $result = PyList_New(len);
  for (i = 0; i < len; i++) {
    PyList_SetItem($result,i,PyString_FromString($1[i]));
    free($1[i]);
  }
}
*/

%{
#include "../src/distorted_sampling.h"
%}

%include "../src/distorted_sampling.h"


