//Para fazer analise: mpip.sourceforge.net

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <mpi.h>

#define MAX_LINHAS 10000
#define MAX_COLUNAS 10000
#define FATOR_ERRO 0.5

double *leArquivo(char *caminho, int linhas, int colunas){
	FILE* arq;
	char ch;
	double *m = malloc(sizeof(double)*linhas*colunas);
	int i, j, k, l;
	arq = fopen(caminho, "r");
	if(arq == NULL){
	    printf("Arquivo não abriu\n");
	}else{
		for(i = 0; i<linhas; i++){
			for (j = 0; j < linhas; j++){
				fscanf(arq,"%lf ",&m[i*colunas+j]);
			}
			fscanf(arq,"%lf \n",&m[i*colunas+j]);
		}
	}

	fclose(arq);
	return m;
}

double *inicializaMatrizU(double *matrizA, int linhas){
	double *matriz = malloc(sizeof(double)*linhas*linhas);
	int i, j;
  for(i=0;i<linhas;i++){
		for(j=0;j<linhas;j++){
      matriz[i*linhas+j] = matrizA[i*linhas+j];
		}
  }
  return matriz;	
}

double *inicializaMatrizL(double *matrizA, int linhas){
	double *matriz = malloc(sizeof(double)*linhas*linhas);
	int i, j;
  for(i=0;i<linhas;i++){
		int j;
		for(j=0;j<linhas;j++){
			if (j==0){
				matriz[i*linhas+j] = matrizA[i*linhas+j];
			}
			if(i==j){
				matriz[i*linhas+j] = 1;
			}
		}
  }
  return matriz;	
}

void imprimeMatriz(double *matriz, int linhas, int colunas){
	int i,j;
  for(i=0;i<linhas;i++){
    for(j=0;j<colunas;j++){
      printf("%.2f ",matriz[i*colunas+j]);
    }
    printf("\n");
  }
	printf("\n");
}


void verificaCorretude(double *matrizL, double *matrizU, double *matriz,int linhas){
	int i,j,k;
	double soma;
	for(i=0;i<linhas;i++){
		for(j=0; j<linhas; j++){
			soma = 0;
			for(k=0; k<linhas; k++){
				soma+= matrizL[i*linhas+k] * matrizU[k*linhas+j];
			}			
			if(soma >= (matriz[i*linhas+j] + FATOR_ERRO) || soma <= (matriz[i*linhas+j] - FATOR_ERRO)){
				printf("Erro na verificacão da linha %d x coluna %d! ",i,j);			
				printf("(esperado: %lf; encontrado %lf)\n",matriz[i*linhas+j],soma);
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

void calculaMatrizLU(double *matrizL, double *matrizU, int linhas, int linhaPivo, int inicio, int fim, int coluna){  
	int i;
		int linhaZerar = inicio;
		for(linhaZerar; linhaZerar < fim; linhaZerar++){
			double coeficiente = (-1)* (matrizU[linhaZerar*linhas+coluna]/matrizU[linhaPivo*linhas+coluna]);
			int j;
			int cont = 0;
			for(j=i; j < linhas; j++){
				matrizU[linhaZerar*linhas+j] = matrizU[linhaPivo*linhas+j] * coeficiente + matrizU[linhaZerar*linhas+j];
				cont++;
			}
			matrizL[linhaZerar*linhas+coluna] = coeficiente * (-1); 
		}
}

int main(int argc, char *argv[]){
	struct struct_info{
		double *matrizU;
		double *matrizL;
		int inicio;
		int fim;
	};


	struct struct_info info;
  int linhas = atoi(argv[2]);
  int colunas = linhas+1;
  double *matriz;
  double *matrizU;
  double *matrizL;

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
		matriz = leArquivo(argv[1],linhas,colunas);
		matrizU = inicializaMatrizU(matriz, linhas);
		matrizL = inicializaMatrizL(matriz, linhas);
		int linhasConsideradas = linhas - 1;
		int divisao;
		int i,j;
		for(i=0;i<(linhas-1);i++){
			divisao = linhasConsideradas/numprocs;
			int qtdProcs = (numprocs>=linhasConsideradas?linhasConsideradas:numprocs);
			linhaPivo = i;
			while(matrizU[linhaPivo*linhas+i] == 0){
				linhaPivo--;
			}

			for(j=1;j<=qtdProcs;j++){
				MPI_Send(&linhaPivo,1,MPI_INT,j,1,MPI_COMM_WORLD);
				MPI_Send(matrizU,(linhas*linhas),MPI_DOUBLE,j,1,MPI_COMM_WORLD);
				MPI_Send(matrizL,linhas*linhas,MPI_DOUBLE,j,1,MPI_COMM_WORLD);
				inicio = fim+1;
				fim = divisao*j;
				if (j == qtdProcs){
					fim = linhas-1;
				}
				MPI_Send(&fim,1,MPI_INT,j,1,MPI_COMM_WORLD);
				MPI_Send(&inicio,1,MPI_INT,j,1,MPI_COMM_WORLD);
			  MPI_Send(&i,1,MPI_INT,j,1,MPI_COMM_WORLD);
			}			
			for(j=1;j<=qtdProcs;j++){
				MPI_Recv(&info,sizeof(info),MPI_PACKED,MPI_ANY_SOURCE,1,MPI_COMM_WORLD, &status);	
				printf("Inicio: %d Fim: %d \n",info.inicio, info.fim);
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
//			matrizU = malloc(sizeof(double)*linhas*linhas);
//			printf("%f\n",matrizU[0]);
			MPI_Recv(matrizU,linhas*linhas,MPI_DOUBLE,mestre,1,MPI_COMM_WORLD, &status);
			printf("%f\n",matrizU[0]);
			MPI_Recv(matrizL,linhas*linhas,MPI_DOUBLE,mestre,1,MPI_COMM_WORLD, &status);
			MPI_Recv(&fim,1,MPI_INT,mestre,1,MPI_COMM_WORLD, &status);
//			printf("Fim: %d\n",fim);
			MPI_Recv(&inicio,1,MPI_INT,mestre,1,MPI_COMM_WORLD, &status);
//			printf("Inicio: %d\n",inicio);
			MPI_Recv(&colunaAtual,1,MPI_INT,mestre,1,MPI_COMM_WORLD, &status);	
//			printf("Coluna: %d\n",colunaAtual);
			calculaMatrizLU(matrizL, matrizU, linhas, linhaPivo, inicio, fim, colunaAtual);
			info.matrizU = matrizU;
			info.matrizL = matrizL;
			info.inicio = inicio;
			info.fim = fim;
			MPI_Send(&info,sizeof(info),MPI_PACKED,mestre,1,MPI_COMM_WORLD);			
		}
  }
  MPI_Finalize();
}
