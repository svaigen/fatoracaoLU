//Para fazer analise: mpip.sourceforge.net

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#define MAX_LINHAS 10000
#define MAX_COLUNAS 10000

float **geraMatrizInicial(int linhas, int colunas){
	float **matriz = malloc(sizeof(float*)*linhas);
	int i, j;
  for(i=0;i<linhas;i++){
    matriz[i] = malloc(sizeof(float)*colunas);
		//for(j=0;j<colunas;j++){
			//matriz[i][j] = (rand() % 10) +1;
		//}
  }
	matriz[0][0] = 5;
	matriz[0][1] = 1;
	matriz[0][2] = -2;
	matriz[0][3] = 10;
	matriz[1][0] = 3;
	matriz[1][1] = -9.4;
	matriz[1][2] = 1.8;
	matriz[1][3] = 22;
	matriz[2][0] = 1;
	matriz[2][1] = 2.2;
	matriz[2][2] = 4.6;
	matriz[2][3] = 10;
  return matriz;
}

float **inicializaMatrizU(float **matrizA, int linhas, int colunas){
	float **matriz = malloc(sizeof(float*)*linhas);
	int i, j;
  for(i=0;i<linhas;i++){
    matriz[i] = malloc(sizeof(float)*(colunas-1));
		for(j=0;j<colunas;j++){
      matriz[i][j] = matrizA[i][j];
		}
  }
  return matriz;	
}

float **inicializaMatrizL(float **matrizA, int linhas, int colunas){
	float **matriz = malloc(sizeof(float*)*linhas);
	int i, j;
  for(i=0;i<linhas;i++){
    matriz[i] = malloc(sizeof(float)*(colunas-1));
		int j;
		for(j=0;j<colunas;j++){
			if (j==0){
				matriz[i][j] = matrizA[i][j];
			}
			if(i==j){
				matriz[i][j] = 1;
			}
		}
  }
  return matriz;	
}

void imprimeMatriz(float **matriz, int linhas, int colunas){
	int i,j;
  for(i=0;i<linhas;i++){
    for(j=0;j<colunas;j++){
      printf("%.2f ",matriz[i][j]);
    }
    printf("\n");
  }
	printf("\n");
}

void **calculaMatrizLU(float **matriz, float **matrizL, float **matrizU, int linhas, int colunas){  
	int i;
	for(i=0; i< colunas-1; i++){
		int linhaZerar = i+1;
		for(linhaZerar; linhaZerar < linhas; linhaZerar++){
			int linhaPivo = linhaZerar-1;
			while(matrizU[linhaPivo][i] == 0){
				linhaPivo--;
			} 
			float coeficiente = (-1)* (matrizU[linhaZerar][i]/matrizU[linhaPivo][i]);
			int j;
			int cont = 0;
			for(j=i; j < colunas; j++){
				matrizU[linhaZerar][j] = matrizU[linhaPivo][j] * coeficiente + matrizU[linhaZerar][j];
				cont++;
			}
			matrizL[linhaZerar][i] = coeficiente * (-1); 
		}
	}
}

float *extraiResultados(float **matriz, int linhas, int colunas){
	int i;
	float *resultados = malloc(sizeof(float)*linhas);
	for(i=0; i< linhas; i++){
		resultados[i] = matriz[i][colunas-1];
	}
	return resultados;
}

void calculaY(float *resultados, float **matrizL, float *vetorY, int linhas){
	vetorY[0] = resultados[0];
	int i,j;
	for(i=1; i < linhas; i++){
		float valor = resultados[i];
		for(j=0; j < i; j++){
			valor -= matrizL[i][j]*vetorY[j];
		}
		vetorY[i] = valor;		
	}	
}

void calculaIncognitas(float *vetorY, float **matrizU, float *incognitas, int linhas, int colunas){
	incognitas[linhas] = vetorY[linhas]/matrizU[linhas][colunas];
	int i,j;
	for(i=linhas-1; i >=0 ; i--){
		float valor = vetorY[i];
		for(j=i+1; j<=colunas; j++){
			valor -= matrizU[i][j] * incognitas[j];
		}
		incognitas[i] = valor / matrizU[i][i];
	} 
}


int main(int argc, char *argv){
  int linhas = 3;
  int colunas = 4;
  float **matriz = geraMatrizInicial(linhas,colunas);
	float **matrizU = inicializaMatrizU(matriz, linhas, colunas-1);
	float **matrizL = inicializaMatrizL(matriz, linhas, colunas-1);

  calculaMatrizLU(matriz, matrizL, matrizU, linhas, colunas-1);

	float *vetorY = malloc(sizeof(double)*linhas);
	float *incognitas = malloc(sizeof(double)*linhas);	
	float *resultados = extraiResultados(matriz,linhas,colunas);
	calculaY(resultados,matrizL,vetorY, linhas);	
	calculaIncognitas(vetorY,matrizU,incognitas, linhas-1, linhas-1);

	printf("Resultado das incognitas: \n\n");
	int i;
	for(i=0; i< linhas; i++){
		printf("%.2f \n",incognitas[i]);
	}
}
