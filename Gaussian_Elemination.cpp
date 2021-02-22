#include <iostream>
#include <fstream> //need to both read and write from/to files.
#include <string>//need to use string to get file names.
#include <stdlib.h> //need to use abs() when calculating condition number and ordering pivots.
using namespace std;
	
void pivotOrganizer(double **A, int n, int i){ //This function checks pivot order for (i+1)th column before every elimination. If needed, it may swap rows.                 
    for (i; i<n; i++)                   //Pivots are ordered respect to their absolute values. Greater absalout value, upper row.
	    for (int k=i+1; k<n; k++)			// Aim to do this order is sending 0 value to the last row.
            if (abs(A[i][i])<abs(A[k][i]))
                for (int j=0;j<=n;j++){//This part swaps rows.
                    float temp=A[i][j];
                    A[i][j]=A[k][j];
                    A[k][j]=temp;}
}

void eliminationF(double **A, int n){//This function does forward Gaussian elimination.
	  for (int i=0; i<n-1; i++){           
		
        	for (int k=i+1; k<n; k++){pivotOrganizer(A, n, i);//This function checks pivot order for (i+1)th column before every reduction. If needed, it may swap rows.
               	 double multiplyingCoefficient=A[k][i]/A[i][i];//Coefficient is calculated by the ratio to pivot
               	 for (int j=0;j<=n;j++)//This loop adds to one row a coefficient multiple of another.
                    A[k][j] -= multiplyingCoefficient*A[i][j]; 
			A[k][i] = 0; // making lower triangle parts exactly zero , eliminating some machine precision errors
				}
    	}
		
		for(int i=0; i<n; i++)//This part checks if matrix is singular, if it is singular, terminates execution
			if(A[i][i] < 0.0000001){// due to machine precision, values which are smaller than 10^-7 are considered as 0
			cout<<"Matrix is singular";
			exit(0);}}



void 	backSubstitution(double **A, int n, double *x){//This function does back substitution and sets x vector with solutions;
	    for (int i=n-1; i>=0; i--){  //Starts from the last row, substitude last x value and goes back to substitude other x values.                      
        	x[i]=A[i][n];                
        	for (int j=i+1;j<n;j++)
            	if (j!=i)            //A[i][i], pivot value, is the desired value                                 
                x[i]=x[i]-A[i][j]*x[j];//All values in row are subtracted except A[i][i], pivot value.
        x[i]=x[i]/A[i][i];            //Now there is only pivot with a coefficient and a value from right hand side. x value is found by division. 
    }
}
    

	
	void solutionPrinter( double *x, int n){//This function creates a file named "Solutions.txt" and writes x values respectively.
	ofstream xFile;
	xFile.open("Solutions.txt");
    for (int i=0; i<n; i++)
	      xFile<< x[i] <<" \n";  
	xFile.close();
	}

 
 
void conditionNumberF(double **A){//This function calculates condition number for 2x2 matrices while detecting if they are singular or not.
	double detA = A[0][0]*A[1][1]-A[0][1]*A[1][0];//Gives determinant of A matrix
	if(detA){//if determinant is not equal to 0, means the matrix not singular. 
	
	//Next line compares sum of absolute values of rows and chooses the greater value, sets to MaxRowSum.
		double MaxRowSum = (abs(A[0][0])+abs(A[0][1])) > (abs(A[1][0])+abs(A[1][1])) ? (abs(A[0][0])+abs(A[0][1])) : (abs(A[1][0])+abs(A[1][1]));
		
	//Next line compares sum of absolute values of columns and chooses the greater value, sets to MaxColumnSum.
		double MaxColumnSum = (abs(A[0][0])+abs(A[1][0])) > (abs(A[0][1])+abs(A[1][1])) ? (abs(A[0][0])+abs(A[1][0])) : (abs(A[0][1])+abs(A[1][1]));
		
		//Condition number can be calculated without taking inverse of the 2x2 matrix.
		//MaxRowSum*MaxColumnSum/detA gives condition number at 1 and infitiy when matrix is 2x2. 
		cout<<"Condition number at 1 and infinity: "<<MaxRowSum*MaxColumnSum/detA;
	
}
	else{
		cout<<"Matrix is singular"<<endl;//When determinant equals to 0.
		exit(0);}
	}



int main(int argc, char *argv[]){
	
	double *b;// b is declared.
	b = new double;
	
	int n;//n is needed during all execution, therefore n is not declared as a pointer
    
    ifstream bFile( argv[2], ios::in );
	if(bFile.is_open()){ // need to check file is open or not.
		while ( bFile>>b[n] )//this loop sets b and get n
		  n++;
	bFile.close();}	
    else
    	cout<<"b file can not open, please re-enter file name"<<endl;

	
	double** A = new double*[n];//A matrix is declared. Later it will be used as augmented matrix.
	for(int i = 0; i < n; i++)
    	A[i] = new double[n+1];//n+1 columns because b will be added to make augmented matrix.
    	
	ifstream AFile( argv[1], ios::in );
	if(AFile.is_open()){ // need to check file is open or not.	
		for(int i=0; i<n; i++){// this loop gets A matrix
			for(int j=0; j<n; j++)
				 AFile>>A[i][j];}
		AFile.close();}		
    else
    	cout<<"A file can not open please re-enter file name"<<endl;	
    	
	//Now A and b matrix are ready to use
	for(int i=0; i<n; i++)//Creating augmented matrix by adding b to end of A
		A[i][n] = b[i]; // A becomes augmented matrix
	//b is not deleted in order to use as x[n] in next steps
	
	if(n==2)
		conditionNumberF(A);//This function prints condition number, and the function may terminate execution if matrix is singular 
		
	eliminationF(A, n);//This function does forward Gaussian elimination, and the function may terminate execution if matrix is singular 
	backSubstitution(A, n, b);//This function does back substitution to get solutions, b is used as x[n] matrix
	solutionPrinter(b, n);//This function writes solutions into "Solutions.txt", b is used as x[n] matrix 

return 0;}


