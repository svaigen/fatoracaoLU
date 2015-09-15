//Para fazer analise: mpip.sourceforge.net

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <time.h>

#define MAX_LINHAS 10000
#define MAX_COLUNAS 10000
#define FATOR_ERRO 0.5

double **leArquivo(char *caminho, int linhas){
	FILE* arq;
	char ch;
	double **m = malloc(sizeof(double*)*linhas);
	int i, j, k, l;
  for(i=0;i<linhas;i++){
    m[i] = malloc(sizeof(double)*(linhas+1));
  }
	arq = fopen(caminho, "r");
	if(arq == NULL){
	    printf("Arquivo não abriu\n");
	}else{
		for(i = 0; i<linhas; i++){
			for (j = 0; j < linhas; j++){
				fscanf(arq,"%lf ",&m[i][j]);
			}
			fscanf(arq,"%lf \n",&m[i][j]);
		}
	}

	fclose(arq);
	return m;
}

double **inicializaMatrizU(double **matrizA, int linhas){
	double **matriz = malloc(sizeof(double*)*linhas);
	int i, j;
  for(i=0;i<linhas;i++){
    matriz[i] = malloc(sizeof(double)*(linhas));
		for(j=0;j<linhas;j++){
      matriz[i][j] = matrizA[i][j];
		}
  }
  return matriz;	
}

double **inicializaMatrizL(double **matrizA, int linhas){
	double **matriz = malloc(sizeof(double*)*linhas);
	int i, j;
  for(i=0;i<linhas;i++){
    matriz[i] = malloc(sizeof(double)*(linhas));
		int j;
		for(j=0;j<linhas;j++){
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

void imprimeMatriz(double **matriz, int linhas, int colunas){
	int i,j;
  for(i=0;i<linhas;i++){
    for(j=0;j<colunas;j++){
      printf("%.2f ",matriz[i][j]);
    }
    printf("\n");
  }
	printf("\n");
}

void calculaMatrizLU(double **matriz, double **matrizL, double **matrizU, int linhas){  
	int i;
	for(i=0; i< linhas; i++){
		int linhaZerar = i+1;
		for(linhaZerar; linhaZerar < linhas; linhaZerar++){
			int linhaPivo = i;
			while(matrizU[linhaPivo][i] == 0){
				linhaPivo--;
			} 
			double coeficiente = (-1)* (matrizU[linhaZerar][i]/matrizU[linhaPivo][i]);
			int j;
			int cont = 0;
			for(j=i; j < linhas; j++){
				matrizU[linhaZerar][j] = matrizU[linhaPivo][j] * coeficiente + matrizU[linhaZerar][j];
				cont++;
			}
			matrizL[linhaZerar][i] = coeficiente * (-1); 
		}
	}
}

double *extraiResultados(double **matriz, int linhas, int colunas){
	int i;
	double *resultados = malloc(sizeof(double)*linhas);
	for(i=0; i< linhas; i++){
		resultados[i] = matriz[i][colunas-1];
	}
	return resultados;
}

void calculaY(double *resultados, double **matrizL, double *vetorY, int linhas){
	vetorY[0] = resultados[0];
	int i,j;
	for(i=1; i < linhas; i++){
		double valor = resultados[i];
		for(j=0; j < i; j++){
			valor -= matrizL[i][j]*vetorY[j];
		}
		vetorY[i] = valor;		
	}	
}

void calculaIncognitas(double *vetorY, double **matrizU, double *incognitas, int linhas, int colunas){
	incognitas[linhas] = vetorY[linhas]/matrizU[linhas][colunas];
	int i,j;
	for(i=linhas-1; i >=0 ; i--){
		double valor = vetorY[i];
		for(j=i+1; j<=colunas; j++){
			valor -= matrizU[i][j] * incognitas[j];
		}
		incognitas[i] = valor / matrizU[i][i];
	} 
}

void verificaCorretude(double **matrizL, double **matrizU, double **matriz,int linhas){
	int i,j,k;
	double soma;
	for(i=0;i<linhas;i++){
		for(j=0; j<linhas; j++){
			soma = 0;
			for(k=0; k<linhas; k++){
				soma+= matrizL[i][k] * matrizU[k][j];
			}			
			if(soma >= (matriz[i][j] + FATOR_ERRO) || soma <= (matriz[i][j] - FATOR_ERRO)){
				printf("Erro na verificacão da linha %d x coluna %d! ",i,j);			
				printf("(esperado: %lf; encontrado %lf)\n",matriz[i][j],soma);
			}
		}
	}
}

void gravaResposta(double *incognitas,char *linhas){
	char caminho[100] = "respostas/\0";
	strcat(caminho,linhas);
	strcat(caminho,".txt");
	FILE* arq;
	arq = fopen(caminho, "w");
	if(arq == NULL){
	    printf("ERRO! Não foi possível criar o arquivo de saída!\n");
	}else{
		int l = atoi(linhas);
		int i;
		for(i = 0; i<l; i++){			
			char valor[20]; 
			sprintf(valor,"%.4f",incognitas[i]);
			fputs(valor,arq);
			fputs("\n",arq);
		}
		printf("Arquivo de saída gerado com sucesso. Os resultados possuem 4 casas decimais. Veja a saída em: %s\n",caminho);
	}

}

int main(int argc, char *argv[]){
	clock_t tempoInicialGeral, tempoFinalGeral, tempoInicialLU, tempoFinalLU;
	tempoInicialGeral = clock();  
	int linhas = atoi(argv[2]);
  int colunas = linhas+1;
  double **matriz = leArquivo(argv[1],linhas);
	double **matrizU = inicializaMatrizU(matriz, linhas);
	double **matrizL = inicializaMatrizL(matriz, linhas);
	tempoInicialLU = clock();
  calculaMatrizLU(matriz, matrizL, matrizU, linhas);
	tempoFinalLU = clock();

	double *vetorY = malloc(sizeof(double)*linhas);
	double *incognitas = malloc(sizeof(double)*linhas);	
	double *resultados = extraiResultados(matriz,linhas,colunas);
	calculaY(resultados,matrizL,vetorY, linhas);	
	calculaIncognitas(vetorY,matrizU,incognitas, linhas-1, linhas-1);
	verificaCorretude(matrizL, matrizU, matriz,linhas);
	gravaResposta(incognitas,argv[2]);
	tempoFinalGeral = clock();
	printf("Tempo de execucao da fatoracao LU: %.8lf segundos \n",(double)(tempoFinalLU - tempoInicialLU)/CLOCKS_PER_SEC);
	printf("Tempo de execucao total: %.8lf segundos \n",(double)(tempoFinalGeral - tempoInicialGeral)/CLOCKS_PER_SEC);
}
