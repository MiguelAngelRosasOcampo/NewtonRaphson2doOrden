#include <stdio.h>
#include <stdlib.h> //Uso de funciones para memoria dinámica.
#include <math.h> //Uso de la función pow.

void newtonRaphson2doOrden(int, double, double*, int*); //Prototipo de función.

int main(void)
{
	system("title Método Newton-Raphson 2do Orden");
	double *Pcoeficientes, terminoInde;
	int *Pexponentes, grado, i;
	
	do 
	{
		setbuf(stdin, NULL); //Limpiar le buffer.
		
		printf("\n\nIngresa el grado del polinomio: ");
		scanf("%d",&grado);
	
	}while(grado < 1 || grado > 50); //Máximo de grado 1 a 50.
	
	Pcoeficientes = (double*) malloc(grado * sizeof(double)); //Reservar memoria dinámica.
	
	Pexponentes = (int*) malloc(grado * sizeof(int));
	
	if(!Pcoeficientes || !Pexponentes) //Verifica que se haya reservado suficiente memoria dinámica.
		exit(-1); //De lo contrario, finalizar el programa.
		
	for(i=0; i<grado; i++) //Ciclo para pedirle al usuario los coeficientes.
	{
		printf("\nIngresa el coeficiente de la variable de grado %d: ",grado-i);
		scanf("%lf",(Pcoeficientes + i));
		*(Pexponentes + i ) = grado - i;
	}
	
	printf("\nAgrega el valor el t%crmino independiente: ",130);
	scanf("%lf",&terminoInde);
	
	newtonRaphson2doOrden(grado, terminoInde, Pcoeficientes, Pexponentes);
	
	free(Pexponentes); //Liberar memoria dinámica.
	free(Pcoeficientes);
	
	getch(); getch(); //Hacer una pausa.
	
	return 0;
}

void newtonRaphson2doOrden(int grado, double terminoInde, double *Pcoeficientes, int *Pexponentes)
{
	double xn = 0, x0, fx = 0, fxD = 0, fxD2 = 0, termino2;
	double termino3, *derivadaCoe, *segundaDerivadaCoe, errorA;
	int i, j = 0, iteraciones;
	
	derivadaCoe = (double*) malloc(grado * sizeof(double)); //Se reserva espacio para los coeficientes.
	segundaDerivadaCoe = (double*) malloc(grado * sizeof(double)); //de las derivadas.
	
	if(!derivadaCoe || !segundaDerivadaCoe)
		exit(-1); //En caso que no se pueda reservar el espacio, salir del programa.
	
	for(i=0; i<grado; i++)//Primera derivada, producto del coeficiente por el exponente de la variable.
		*(derivadaCoe + i) = *(Pcoeficientes + i) * *(Pexponentes + i);
	
	termino2 = *(Pcoeficientes + (i-1)); //Se guarda el término independiente de la primera derivada.
	
	for(i=0; i<grado; i++) //Segunda derivada
		*(segundaDerivadaCoe + i) = *(derivadaCoe + i) * ( *(Pexponentes + i) - 1);
		
	termino3 = *(derivadaCoe + (i-2)); //Se guarda el término independiente de la segunda derivada.
		
	do //Bucle para que usuario ingrese cantidad de iteraciones positivas.
	{
		setbuf(stdin, NULL);
		
		printf("\nIngresa la cantidad de iteraciones: ");
		scanf("%d",&iteraciones);
	
	}while(iteraciones < 1);

	printf("\nIngresa el valor inicial de x0: ");
	scanf("%lf",&x0);
	
	printf("\n\n\t\t\t\tM%ctodo Newton-Raphson 2do Orden\n\n\t\t\t\t%c(x)   = ",130,159); //Mostrar F(x)
	
	for(i=0; i<grado; i++) //Bucle para mostrarle la función al usuario.
	{
		if( !*(Pcoeficientes + i)  ) //Si el coeficiente es 0, NO mostrar.
			continue;
			
		else printf("%+5gx^%d ",*(Pcoeficientes+i), *(Pexponentes+i));
	}

	printf("%+5g",terminoInde); //Mostrar el término independiente.
	
	printf("\n\n\t\t\t\t%c\'(x)  = ",159); //Mostrar F'(x).

	for(i=0; i<grado-1; i++) //Mostrar la derivada de la función.
	{
		if( !*(derivadaCoe + i) ) //Si el coeficiente es 0, NO mostrar.
			continue;
			
		else printf("%+5gx^%d ",*(derivadaCoe + i), *(Pexponentes + i) - 1);
	}
		
	printf("%+5g",termino2); //Mostrar el término independiente de la derivada.

	printf("\n\n\t\t\t\t%c\'\'(x) = ",159); //Mostrar F''(x).

	for(i=0; i<grado-2; i++)
	{
		if( !*(segundaDerivadaCoe + i) ) //Si el coeficiente es 0, NO mostrar.
			continue;
			
		else printf("%+5gx^%d ",*(segundaDerivadaCoe + i), *(Pexponentes + i) - 2);
	}

	printf("%+5g",termino3); //Mostrar el término independiente de la segunda derivada.

	printf("\n\n\n I \t\tXn \t\tXn+1 \t\t%c(x) \t\t%c\'(x) \t\t%c\'\'(x) \t\teA",159,159,159);
	
	while(j++ != iteraciones)
	{			
		fx = 0;		fxD = 0;	fxD2 = 0; //Las variables se inicializan en 0 en cada iteración.

		for(i=0; i<grado; i++) //Obtener el valor de la evaluación de f(x).
			fx += *(Pcoeficientes+i) * pow(x0, *(Pexponentes+i));
		
		fx += terminoInde; //Sumar el término independiente.
		
		for(i=0; i<grado-1; i++) //Obtener el valor de la evaluación de f'(x).
			fxD += *(derivadaCoe + i) * pow(x0, *(Pexponentes + i ) - 1);
			
		fxD += termino2; //Sumar el término independiente.
		
		for(i=0; i<grado-2; i++) //Obtener el valor de la evaluación de f''(x).
			fxD2 += *(segundaDerivadaCoe + i) * pow(x0, *(Pexponentes + i) - 2);
		
		fxD2 += termino3; //Sumar el término independiente.		

		xn = x0 - (fx*fxD)/(pow(fxD,2) - fx*fxD2); //Calcular Xn+1
		errorA = fabs( (xn - x0) / xn) * 100; //Obtener el error.
				
		printf("\n\n %d \t\t%lf \t%lf \t%lf \t%lf \t%lf \t%lf %%",j,x0,xn,fx,fxD,fxD2,errorA);		

		x0 = xn;		
	}
	
	free(derivadaCoe); //Liberar memoria dinámica.
	free(segundaDerivadaCoe);
}
