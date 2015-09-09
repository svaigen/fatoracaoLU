//Para fazer analise: mpip.sourceforge.net

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <mpi.h>

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

void calculaMatrizLU(double **matrizL, double **matrizU, int linhas, int linhaPivo, int inicio, int fim, int coluna){  
	int i;
		int linhaZerar = inicio;
		for(linhaZerar; linhaZerar < fim; linhaZerar++){
			double coeficiente = (-1)* (matrizU[linhaZerar][coluna]/matrizU[linhaPivo][coluna]);
			int j;
			int cont = 0;
			for(j=i; j < linhas; j++){
				matrizU[linhaZerar][j] = matrizU[linhaPivo][j] * coeficiente + matrizU[linhaZerar][j];
				cont++;
			}
			matrizL[linhaZerar][coluna] = coeficiente * (-1); 
		}
}

int main(int argc, char *argv[]){
  int linhas = atoi(argv[2]);
  int colunas = linhas+1;
  double **matriz;
  double **matrizU;
  double **matrizL;

  int myid, numprocs;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	numprocs--;
	int mestre = 0;
	int linhaPivo;
	int inicio, fim = 0;

  if(myid == 0){
		matriz = leArquivo(argv[1],linhas);
		matrizU = inicializaMatrizU(matriz, linhas);
		matrizL = inicializaMatrizL(matriz, linhas);
		int linhasConsideradas = linhas - 1;
		int divisao;
		int i,j;
		for(i=0;i<(linhas-1);i++){
			divisao = linhasConsideradas/numprocs;
			int qtdProcs = (numprocs>=linhasConsideradas?numprocs:linhasConsideradas);
			linhaPivo = i;
			while(matrizU[linhaPivo][i] == 0){
				linhaPivo--;
			}

			for(j=1;j<=qtdProcs;j++){
				MPI_Send(&linhaPivo,1,MPI_INT,j,1,MPI_COMM_WORLD);
				MPI_Send(&matrizU,linhas*linhas,MPI_DOUBLE,j,1,MPI_COMM_WORLD);
				MPI_Send(&matrizL,linhas*linhas,MPI_DOUBLE,j,1,MPI_COMM_WORLD);
				inicio = fim+1;
				fim = divisao*j;
				if (j == qtdProcs){
					fim = linhas;
				}
				MPI_Send(&inicio,1,MPI_INT,j,1,MPI_COMM_WORLD);
				MPI_Send(&fim,1,MPI_INT,j,1,MPI_COMM_WORLD);
			  MPI_Send(&i,1,MPI_INT,j,1,MPI_COMM_WORLD);
			}
			
			linhasConsideradas--;
		}

  }else{
		int colunaAtual;
		while(1){
			MPI_Recv(&linhaPivo,1,MPI_INT,mestre,1,MPI_COMM_WORLD, &status);
			if (linhaPivo == -1){
				break;
			}
			MPI_Recv(&matrizU,1,MPI_DOUBLE,mestre,1,MPI_COMM_WORLD, &status);
			MPI_Recv(&matrizL,1,MPI_DOUBLE,mestre,1,MPI_COMM_WORLD, &status);
			MPI_Recv(&fim,1,MPI_INT,mestre,1,MPI_COMM_WORLD, &status);
			MPI_Recv(&inicio,1,MPI_INT,mestre,1,MPI_COMM_WORLD, &status);
			MPI_Recv(&colunaAtual,1,MPI_INT,mestre,1,MPI_COMM_WORLD, &status);	

			calculaMatrizLU(matrizL, matrizU, linhas, linhaPivo, inicio, fim, coluna);
			
			MPI_Send(&matrizU,linhas*linhas,MPI_DOUBLE,mestre,1,MPI_COMM_WORLD);
			MPI_Send(&matrizL,linhas*linhas,MPI_DOUBLE,mestre,1,MPI_COMM_WORLD);
			MPI_Send(&inicio,1,MPI_INT,mestre,1,MPI_COMM_WORLD);
		  MPI_Send(&fim,1,MPI_INT,mestre,1,MPI_COMM_WORLD);
		}

  }
  MPI_Finalize();
}
