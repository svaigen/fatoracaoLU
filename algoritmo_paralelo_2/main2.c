//Para fazer analise: mpip.sourceforge.net

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <mpi.h>
#include <time.h>

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
	int soma =0;
  for(i=0;i<linhas;i++){
		for(j=0;j<linhas;j++){
      matriz[i*linhas+j] = matrizA[i*linhas+j+soma];
		}
		soma++;
  }
  return matriz;	
}

double *inicializaMatrizL(double *matrizA, int linhas){
	double *matriz = malloc(sizeof(double)*linhas*linhas);
	int i, j;
  for(i=0;i<linhas;i++){
		int j;
		for(j=0;j<linhas;j++){
			if(i==j){
				matriz[i*linhas+j] = 1;
			}else{
				matriz[i*linhas+j] = 0;
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
	int correcao = 0;
	for(i=0;i<linhas;i++){
		for(j=0; j<linhas; j++){
			soma = 0;
			for(k=0; k<linhas; k++){
				soma+= matrizL[i*linhas+k] * matrizU[k*linhas+j];
			}			
			if(soma >= (matriz[i*linhas+j+correcao] + FATOR_ERRO) || soma <= (matriz[i*linhas+j+correcao] - FATOR_ERRO)){
				printf("Erro na verificacão da linha %d x coluna %d! ",i,j);			
				printf("(esperado: %lf; encontrado %lf)\n",matriz[i*linhas+j+correcao],soma);
			}
		}
		correcao++;
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

double *extraiResultados(double *matriz, int linhas, int colunas){
	int i;
	double *resultados = malloc(sizeof(double)*linhas);
	for(i=0; i< linhas; i++){
		resultados[i] = matriz[i * colunas +colunas-1];
	}
	return resultados;
}

void calculaY(double *resultados, double *matrizL, double *vetorY, int linhas){
	vetorY[0] = resultados[0];
	int i,j;
	for(i=1; i < linhas; i++){
		double valor = resultados[i];
		for(j=0; j < i; j++){
			valor -= matrizL[i*linhas+j]*vetorY[j];
		}
		vetorY[i] = valor;		
	}	
}

void calculaIncognitas(double *vetorY, double *matrizU, double *incognitas, int linhas, int colunas){
	incognitas[linhas] = vetorY[linhas]/matrizU[linhas*colunas-1];
	int i,j;
	for(i=linhas-1; i >=0 ; i--){
		double valor = vetorY[i];
		for(j=i+1; j<colunas; j++){
			valor -= matrizU[i*colunas+j] * incognitas[j];
		}
		incognitas[i] = valor / matrizU[i*linhas+i];
	} 
}


void calculaMatrizLU(double *matrizL, double *matrizU, int linhas, int linhaPivo, int inicio, int fim, int coluna){  
		int linhaZerar = inicio;
		for(linhaZerar; linhaZerar < fim+1; linhaZerar++){
			double coeficiente = (-1)* (matrizU[linhaZerar*linhas+coluna]/matrizU[linhaPivo*linhas+coluna]);
			int j;
			for(j=0; j < linhas; j++){
				matrizU[linhaZerar*linhas+j] = matrizU[linhaPivo*linhas+j] * coeficiente + matrizU[linhaZerar*linhas+j];
			}
			matrizL[linhaZerar*linhas+coluna] = coeficiente * (-1); 
		}
}

void criaNovaMatrizLU(double *matrizL, double *matrizU, int linhas, int inicio, int fim, int coluna, double *novaMatrizU, double *novaMatrizL){  
	int i,j;
		for(i=inicio;i<fim+1;i++){
			for(j=0;j<linhas;j++){
				matrizU[i*linhas+j] = novaMatrizU[i*linhas+j];
			}
			matrizL[i*linhas+coluna] = novaMatrizL[i*linhas+coluna];	
		}
}

void copiaMatrizLU(double *matrizL, double *matrizU, int linhas, int inicio, int fim, int coluna, double *infoMatrizU, double *infoMatrizL){  
	int i,j;
		for(i=inicio;i<fim+1;i++){
			for(j=0;j<linhas;j++){
				infoMatrizU[i*linhas+j] = matrizU[i*linhas+j];
			}
			infoMatrizL[i*linhas+coluna] = matrizL[i*linhas+coluna];	
		}
}

int main(int argc, char *argv[]){
  clock_t tempoInicialGeral, tempoFinalGeral, tempoInicialLU, tempoFinalLU;

	tempoInicialGeral = clock();

	int i;
  int linhas = atoi(argv[2]);
  int colunas = linhas+1;
  double *matriz;
  double *matrizU;
  double *matrizL;

	matrizU = malloc(sizeof(double)*linhas*linhas-1);
	matrizL = malloc(sizeof(double)*linhas*linhas-1);
  int myid, numprocs;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	int escravos = numprocs-1;
/*operacoes referentes ao broadcast*/
	int destinos[escravos];
	for(i=1; i< numprocs; i++){
		destinos[i-1] = i;
	}
  MPI_Group grupoGlobal;
  MPI_Comm_group(MPI_COMM_WORLD, &grupoGlobal);
	MPI_Group grupoEscravos;
  MPI_Group_incl(grupoGlobal, (escravos), destinos, &grupoEscravos);
	MPI_Comm commEscravos;
  MPI_Comm_create(MPI_COMM_WORLD, grupoEscravos, &commEscravos);	

	int totalElementos = linhas*linhas;
	int mestre = 0;
	int linhaPivo;
	int inicio, fim = 0;

	int qtdProcs;

  if(myid == mestre){
		matriz = leArquivo(argv[1],linhas,colunas);
		matrizU = inicializaMatrizU(matriz, linhas);
		matrizL = inicializaMatrizL(matriz, linhas);
		//novaU = inicializaMatrizU(matriz, linhas);
		//novaL = inicializaMatrizU(matriz, linhas);
    MPI_Bcast(matrizU,totalElementos,MPI_DOUBLE,mestre,MPI_COMM_WORLD);		
    MPI_Bcast(matrizL,totalElementos,MPI_DOUBLE,mestre,MPI_COMM_WORLD);		

		int linhasConsideradas = linhas - 1;
		int divisao;
		int j;
		tempoInicialLU = clock();
		for(i=0;i<(linhas-1);i++){
			divisao = linhasConsideradas/escravos;
			qtdProcs = (escravos>=linhasConsideradas?linhasConsideradas:escravos);
			linhaPivo = i;
			while(matrizU[linhaPivo*linhas+i] == 0){
				linhaPivo--;
			}
			if (divisao == 0){	
				inicio = i+1;
			}

//		    MPI_Bcast(&linhaPivo,1,MPI_INT,mestre,MPI_COMM_WORLD);		
//				MPI_Bcast(&i,1,MPI_INT,mestre,MPI_COMM_WORLD);		
//		    MPI_Bcast(&qtdProcs,1,MPI_INT,mestre,MPI_COMM_WORLD);		

			for(j=1;j<=qtdProcs;j++){
				MPI_Send(&linhaPivo,1,MPI_INT,j,1,MPI_COMM_WORLD);				
				MPI_Send(&i,1,MPI_INT,j,2,MPI_COMM_WORLD);		
		    MPI_Send(&qtdProcs,1,MPI_INT,j,3,MPI_COMM_WORLD);		

				if (divisao == 0){
					fim = inicio;	
				}else{
				inicio = (j-1) * divisao + i;
					inicio++;
				fim = j*divisao +i;
				}
				if (j == qtdProcs){
					fim = linhas-1;
				}
				MPI_Send(&fim,1,MPI_INT,j,4,MPI_COMM_WORLD);
				MPI_Send(&inicio,1,MPI_INT,j,5,MPI_COMM_WORLD);
				inicio++;
			}			
//			MPI_Recv(matrizU,totalElementos,MPI_DOUBLE,MPI_ANY_SOURCE,1,MPI_COMM_WORLD, &status);
//			MPI_Recv(matrizL,totalElementos,MPI_DOUBLE,MPI_ANY_SOURCE,1,MPI_COMM_WORLD, &status);
			linhasConsideradas--;
		}
		tempoFinalLU=clock();
		MPI_Recv(matrizU,totalElementos,MPI_DOUBLE,MPI_ANY_SOURCE,1,MPI_COMM_WORLD, &status);
		MPI_Recv(matrizL,totalElementos,MPI_DOUBLE,MPI_ANY_SOURCE,1,MPI_COMM_WORLD, &status);
		//Encerra outros processos
		linhaPivo = -1;
    //MPI_Bcast(&linhaPivo,1,MPI_INT,mestre,MPI_COMM_WORLD);		
				
		for(i=1;i<=escravos;i++){
			MPI_Send(&linhaPivo,1,MPI_INT,i,1,MPI_COMM_WORLD);
		}


		double *vetorY = malloc(sizeof(double)*linhas);
		double *incognitas = malloc(sizeof(double)*linhas);	
		double *resultados = extraiResultados(matriz,linhas,colunas);
		calculaY(resultados,matrizL,vetorY, linhas);	
		calculaIncognitas(vetorY,matrizU,incognitas, linhas, linhas);
		verificaCorretude(matrizL, matrizU, matriz,linhas);
		gravaResposta(incognitas,argv[2]);
		tempoFinalGeral = clock();
		printf("Tempo de execucao da fatoracao LU: %.8lf segundos \n",(double)(tempoFinalLU - tempoInicialLU)/CLOCKS_PER_SEC);
		printf("Tempo de execucao total: %.8lf segundos \n",(double)(tempoFinalGeral - tempoInicialGeral)/CLOCKS_PER_SEC);
  }else{
		MPI_Bcast(matrizU,totalElementos,MPI_DOUBLE,mestre,MPI_COMM_WORLD);		
    MPI_Bcast(matrizL,totalElementos,MPI_DOUBLE,mestre,MPI_COMM_WORLD);
		int colunaAtual;
		while(1){
	    MPI_Recv(&linhaPivo,1,MPI_INT,mestre,1,MPI_COMM_WORLD,&status);							
			if (linhaPivo == -1){
				break;
			}
			MPI_Recv(&colunaAtual,1,MPI_INT,mestre,2,MPI_COMM_WORLD,&status);		
	    MPI_Recv(&qtdProcs,1,MPI_INT,mestre,3,MPI_COMM_WORLD,&status);
			
			MPI_Recv(&fim,1,MPI_INT,mestre,4,MPI_COMM_WORLD, &status);
			MPI_Recv(&inicio,1,MPI_INT,mestre,5,MPI_COMM_WORLD, &status);
			//MPI_Recv(&matrizU[inicio*linhas],((fim-inicio+1)*linhas),MPI_DOUBLE,mestre,1,MPI_COMM_WORLD, &status);
			//MPI_Recv(&matrizL[inicio*linhas],((fim-inicio+1)*linhas),MPI_DOUBLE,mestre,1,MPI_COMM_WORLD, &status);


			calculaMatrizLU(matrizL, matrizU, linhas, linhaPivo, inicio, fim, colunaAtual);

			//momento de envio para todos			
			int i;
			for (i = 0; i < qtdProcs; i++){
	      int posInicial = inicio*linhas;
  	    int qtdeEnviados = ((fim-inicio+1)*linhas);
				MPI_Bcast(&posInicial,1,MPI_INT,i,commEscravos);
        MPI_Bcast(&qtdeEnviados,1,MPI_INT,i,commEscravos);
        MPI_Bcast(&matrizU[posInicial],qtdeEnviados,MPI_DOUBLE,i,commEscravos);
        MPI_Bcast(&matrizL[posInicial],qtdeEnviados,MPI_DOUBLE,i,commEscravos);
      }

			if (myid == 1 && qtdProcs == 1){
    	  MPI_Send(matrizU, totalElementos, MPI_DOUBLE, mestre, 1, MPI_COMM_WORLD);
    	  MPI_Send(matrizL, totalElementos, MPI_DOUBLE, mestre, 1, MPI_COMM_WORLD);
    	}

			//MPI_Send(&fim,1,MPI_INT,mestre,1,MPI_COMM_WORLD);
			//MPI_Send(&inicio,1,MPI_INT,mestre,1,MPI_COMM_WORLD);			
			//MPI_Send(&matrizU[inicio*linhas],((fim-inicio+1)*linhas),MPI_DOUBLE,mestre,1,MPI_COMM_WORLD);
			//MPI_Send(&matrizL[inicio*linhas],((fim-inicio+1)*linhas),MPI_DOUBLE,mestre,1,MPI_COMM_WORLD);			
		}
  }
  MPI_Finalize();
}
