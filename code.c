#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "eispack.h"
#define  MAX 100

//160202093 Muhammed Ali Dilekci II. Ogretim

void ozvektorbul()
{
    int boyut;
    FILE *dosya;
    dosya = fopen("ozvektorler.txt","a");


    printf("nxn boyutlu matris icin n giriniz:"); scanf("%d",&boyut);
    int i;
    double* a = ( double * ) malloc (boyut * boyut * sizeof ( double ) );
    double* w = ( double * ) malloc (boyut * sizeof ( double ) );
    double* x = ( double * ) malloc (boyut * boyut * sizeof ( double ) );
    for(i = 0; i < boyut*boyut;i++)
    {
        double v;
        printf(" Matris[%d][%d] : ", i/boyut + 1, i%boyut + 1);
        scanf("%lf", &v);
        a[i] = v;

    }
	int matz = 1;
	int ierr = rs ( boyut, a, w, matz, x );

    fprintf(dosya,"\nOzdegerler : \n");
    printf("Ozdegerler : \n");
    for(i=0;i<boyut;i++)
    {
        fprintf(dosya,"%.2f",w[i]);
        printf("%.2f ",w[i]);
    }

    printf("\nOzvektorler : \n");
    for ( i = 0; i < boyut; i++ )
    {
        fprintf(dosya,"%.2f",x[i]);
        printf("%.2f ",x[i]);
    }
}

void nilpotent()
{
FILE *dosya;
    dosya=fopen("nilpotent.txt","a");
int m=0,n=0,i=0,j=0,toplam=0;

    printf(" nxn matris icin n giriniz : ");    scanf("%d",&n);

int carpim[n][n],matris[n][n];
    srand(time(NULL));

    for(i=0; i<n; i++)                      //Matrise rastgele deger atama
    {
        for(j=0; j<n; j++)
        {
            matris[i][j]=(rand()%30)-15;
        }
    }

    printf("Rastgele matrisiniz\n");        //Rastgele matrisi ekrana yazdirma
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            printf(" %d ",matris[i][j]);
        }
        printf("\n");
    }

    fprintf(dosya,"\nRastgele  matris : \n");//Rastgele matrisi dosyaya yazma
    for(i=0; i<n; i++)
    {
        fprintf(dosya,"|");
        for(j=0; j<n; j++)
            fprintf(dosya,"%d ",matris[i][j]);
        fprintf(dosya,"|\n");
    }
    printf("\n");

    for( i=0; i<n; i++)                     //Matrisleri carpma
    {
        for( j=0; j<n; j++)
        {
            for( m=0; m<n; m++)
            {
                toplam += (matris[i][m] * matris[m][j]);
            }
            carpim[i][j] = toplam;
            toplam = 0;
        }
        printf("\n");
    }
    printf("\n");

    fprintf(dosya,"\nOlusturulan  matris : \n");//Ic carpilmis matrisi dosyaya yazma
    for(i=0; i<n; i++)
    {
        fprintf(dosya,"|");
        for(j=0; j<n; j++)
            fprintf(dosya,"%d ",carpim[i][j]);
            fprintf(dosya,"|\n");
    }

int test=0;                                     // Toplam==0 kontrolu
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            if (carpim[i][j]!=0) test=1;
        }
    }

    if (test)
    {
        fprintf(dosya,"\nNilpotent degildir.\n");
        printf("Nilpotent degildir.");
    }
    if(!test)
    {
        fprintf(dosya,"\nNilpotenttir.\n");
        printf("Nilpotenttir.");
    }
}

schur(double w[], double a[],int n)       //w => ozdegerler   a => elemanlar
{
double lmnKare=0;
double ozKare=0;
FILE *dosya;
dosya=fopen("Schur.txt","a");
int i;
    for(i=0; i<n*n; i++)      //Elemanlarin kareleri toplamini bulur.
        lmnKare+=a[i]*a[i];

    for(i=0; i<n; i++)        //ozDegerlerin kareleri toplamini bulur.
        ozKare +=w[i]*w[i];

    fprintf(dosya,"\n Elemanlar : \n|");    //Elemanlari   dosyaya yazma
    for(i=0; i<n*n; i++)
    {
        if(i%n==0 && i!=0)fprintf(dosya,"|\n|");
        fprintf(dosya,"%.2f ",a[i]);

    }
    fprintf(dosya,"|");

    fprintf(dosya,"\n Ozdegerler : \n");    //ozDegerleri  dosyaya yazma
    for(i=0; i<n; i++)
    {
        fprintf(dosya,"%.2f ",w[i]);
    }

    fprintf(dosya,"\n Ozdegerlerin kareleri toplami : %.2f\n",ozKare);    //Schur karsilastirmalarini dosyaya yazar
    fprintf(dosya,"\n Elemanlarin  kareleri toplami : %.2f\n",lmnKare);
    if(ozKare<lmnKare) fprintf(dosya,"\nSchur esitsizligine uyar.\n");
    else fprintf(dosya,"\nSchur esitsizligine uymaz.\n");

    printf("\nozdegerlerin kareleri toplami: %.2f\n",ozKare);    //Schur karsilastirmalarini ekrana yazar
    printf("\nElemanlarin kareleri toplami: %.2f\n",lmnKare);
    if(ozKare<lmnKare) printf("\nSchur esitligine uyar.\n");
    else printf("\nSchur esitsizligine uymaz.\n");
}

int ozdeger(int disli)      // disli :Ozdeger isleminde 1 schur isleminde 2 olur
{
    FILE *dosya;
    dosya=fopen("ozdegerler.txt","a");
int n;
    printf("nxn matris icin n giriniz : ");     // Matris boyutunu alma
    scanf("%d",&n);
int i;
double delta=0;
double matris[n][n];
double ozDegerArr[n];

    for(i=0; i<n; i++)ozDegerArr[i]=0;          //Ozdegerleri tutacak matrisin butun elemanlarini 0 a esitler

    if(n==2)                                    //2x2 matrislerin ozdegerini bulmak
    {
double a,b,c,k,l,m,n;                           //Matris Elemanlarini alma
        printf("1.satir 1.sutun elemanini giriniz.");
        scanf("%lf",&matris[0][0]);
        printf("1.satir 2.sutun elemanini giriniz.");
        scanf("%lf",&matris[0][1]);
        printf("2.satir 1.sutun elemanini giriniz.");
        scanf("%lf",&matris[1][0]);
        printf("2.satir 2.sutun elemanini giriniz.");
        scanf("%lf",&matris[1][1]);
        k=matris[0][0];
        l=matris[0][1];
        m=matris[1][0];
        n=matris[1][1];

        fprintf(dosya,"\nGirilen matris :\n|%.2f ,%.2f| \n|%.2f %.2f| ",matris[0][0],matris[0][1],matris[1][0],matris[1][1]);//Girilen matrisi yazdirma

        a=1;                            // ax^2+bx+c=0 icin a,b,c,delta hesaplama
        b=(-k-n);
        c=(k*n)-(l*m);
        delta=(double)(b*b)-(4*a*c);

        if(delta<0)             //Denklem kokleri= ozdegerler
        {
            printf("Reel ozdegerler bulunmamaktadir.");
            return 0;
        }

        if(delta ==0)
        {
            ozDegerArr[0]= (b)/(-2);
        }
        if(delta>0)
        {
            ozDegerArr[0]= ((-1*b)+ pow(delta,0.5))/2;
            ozDegerArr[1]= ((-1*b)- pow(delta,0.5))/2;
        }


        if(disli==1)            //Ozdeger islemi ise     ozdegerleri ekrana ve ozdeger.txt ye yazdir
            printf("\nozdegerler : %f  %f ",ozDegerArr[0],ozDegerArr[1]);

        fprintf(dosya,"\nozdegerler :%f ,%f",ozDegerArr[0],ozDegerArr[1]);

        if(disli==2)            //Schur islemi ise
        {
            schur(ozDegerArr,matris,2);
        }
    }

    if(n>2)                 //2x2 den buyukse
    {


int matz = 1;

double* a = ( double * ) malloc (n * n * sizeof ( double ) );
double* w = ( double * ) malloc (n * sizeof ( double ) );
double* x = ( double * ) malloc (n * n * sizeof ( double ) );

int i;
        for( i = 0; i < n*n; i++)   //Matris elemanlarini almak
        {
            double v;
            printf(" Matris[%d][%d] : ", i/n + 1, i%n + 1);
            scanf("%lf", &v);
            a[i] = v;
        }

int ierr = rs ( n, a, w, matz, x );//Ozdeger hesaplayip w[] dizisine atar. matris elemanlarini a[] dizisine atar

        if ( ierr != 0 )                   // Hata kontrolu
        {
            printf ( "\nHata :(\n" );
            return;
        }
        if(disli!=2)                       // Ozdeger islemiyse(ozdegerleri yazdir)
        {
            fprintf(dosya,"\nGirilen Matris :\n|"); //Girilen matrisi ekrana ve dosyaya yazdirma
            printf("\nGirilen Matris :\n|");
            for(i=0; i<n*n; i++)
            {
                fprintf(dosya,"%.2f ",a[i]);
                printf("%.2f ",a[i]);
                if( ((i+1)%n==0) && ((i+1)!=n*n) && (i!=0) )
                {
                    fprintf(dosya,"|\n|");
                    printf("|\n|");
                }
                if ( i==(n*n)-1 )
                {
                    fprintf(dosya,"|");
                    printf("|");
                }
            }

            fprintf(dosya,"\nOzdegerler :");        //Ozdegerleri ekrana ve dosyaya yazdirma
            printf("\nOzdegerler:\n");
            for (i = 0; i < n; i++ )
            {
                printf (" => %f\n", w[i]);
                fprintf(dosya,"%f ",w[i]);
            }
        }

        if (disli==2)                      // Schur islemiyse (ozDegerleri ve elemanlari schur a yolla)
        {
double elemanlar[n];
double ozdegerler[n];
            for (i=0; i<n; i++)   ozdegerler[i]=w[i];
            for (i=0; i<n*n; i++)elemanlar[i]=a[i];
            schur(ozdegerler,elemanlar,n);
        }
    }

    return 0;
}

int cark(int kac)
{
int i,counter;
srand(time(NULL));

    for(i=1; i<=kac; i++)       // 0- 241 arasi random deger alir
    {
        counter=rand()%242;
        printf("%d   ",counter);
        if (i%10==0) printf("\n");
    }
    printf("\n\n %d carktan cikan son sayidir.",counter);
    return counter;
}


int main()
{
int kac,disli;

    printf("Carkin kac kez donmesini istediginizi giriniz : ");
    scanf("%d",&kac);
    disli=cark(kac);
    switch(disli)
    {

    case 1:
        printf("Ozdeger islemindesiniz.\n\n");
        ozdeger(1);
        break;
    case 2:
        printf("Schur islemindesiniz.\n\n");
        ozdeger(2);
        break;
    case 3:
        printf("Ozvektor islemindesiniz.\n\n");
        ozvektorbul();
        break;
    case 4:
        printf("Nilpotent islemindesiniz.\n\n");
        nilpotent();
        break;
    }


    return 0;
}
