/* 
27.3.2020: Corrected error that could lead to the generation of non-isomorphic
graphs. No graphs were lost due to this error, but isomorphic copies could be generated.
The error was in a part that was only used if the degree sequence prescibed at least 7 
vertices of degree one and sufficiently many vertices of degree at least 3 that the
graphs were no trees.

In addition, this version (the first modified one since about 25 years...) also contains
modifications by Brendan McKay who changed the old-style function declarations
to new ones, so that gcc version 9 complains a bit less...
*/

/*****MULTIGRAPH.C*********/
/* Unterschied zu minimulti.c: Es wird daraufhingearbeitet verschiedene
   Restriktionen, wie z.B. vorgegebene Taillenweite, ... mit einzubeziehen */
/*

Generate connected multigraphs.

The first parameters are the degree sequence, starting with the number of
   vertices of degree 1, then 2, etc

Options are "o" -- output the graphs in multicode -- default stdout,
   but to a file with default name if you use option "f".
For counting only, omit both "o" and "f".

"tx" -- generate graphs with girth at least x (so t3 for simple graphs)

"b" for bipartite graphs

"m x y" with 0<=x<y for the usual modulo splitting

So "multigraph 0 3 4 11 t6 b o" would generate all graphs with 3 vertices
   of degree 2, 4 of degree 3, 11 of degree 4 that are bipartite, have
   girth at least 6. The graphs would be written in multicode to stdout.
*/

# include <stdio.h>
# include <stdlib.h>
//# include <malloc.h>
# include <time.h>
# include <string.h>
# include <sys/times.h>
# include "multigraph.h"

/*******DEFINITIONEN***************************************/

/*******Should*be*changed:********************************/

# define knoten 63
  /* maximal value 64 -- at least one larger than the used vertexnumber*/ 

# define baumgrenze 40

#define listenlaenge 10000000

/******Should*not*be*changed:*****************************/

#define Sicherungsintervall 7200

                /*At least "Sicherungsintervall" seconds between two savings*/

/* wird momentan nicht benutzt -- muss auf signale wie in minibaum umgestellt werden */

/***********May*not*be*changed:*****************************/


#define unbelegt 98

#define leer 99

#define nil 0

#if !defined(CLK_TCK) && !defined(_SC_CLK_TCK)
#include <unistd.h>
#endif
#if !defined(CLK_TCK) && defined(_SC_CLK_TCK)
#define CLK_TCK sysconf(_SC_CLK_TCK)
#endif
#if !defined(CLK_TCK) && defined(CLOCKS_PER_SEC)
#define CLK_TCK CLOCKS_PER_SEC
#endif

#ifdef CLK_TCK
#define time_factor CLK_TCK
#endif


typedef unsigned char GRAPH[knoten+1][knoten+1];


typedef unsigned char BILD[knoten+1];
typedef unsigned char ADJAZENZ[knoten+1];
typedef unsigned char BILDELEMENT;


typedef struct zweig {
             signed char blatt;            /* 1, wenn es ein blatt ist, 0 sonst */
             BILD bild;
             unsigned char *tot;    /* Ein array mit Totwerten. tot[i] 
				       gibt an, ob es noch Sinn hat, in die 
				       i-te Richtung weiterzusuchen, oder ob da
                                       keine Moeglichkeit besteht */
             struct zweig ** successors;   /* Die Nachfolger */
             int anzahlnachfolger;
             int anzahltot;         /* Die Anzahl der "toten" Nachfolger */
             struct zweig *vorgaenger;
             unsigned char rootnumber; /* Welche nummer hat die ROOT dieses
                                          Baumes -- oder auch: Was ist das
                                          Urbild von 1 */
             int zweignummer;       /* Der wievielte nachfolger des vorigen
                                       Zweiges ist dieser */
             unsigned char nextbild;/* die naechste zu vergebende nummer */
             unsigned char wo;      /* welchen "wahren" Namen hat der naechste
                                       Knoten, an dem noch wahlmoeglichkeiten
                                       fuer bild bestehen, oder der noch nicht
                                       fertig ist */
             unsigned char tiefe;   /* in welcher tiefe wurde dieser knoten
                                       eingefuegt (d.h.: Wie gross ist die
                                       anzahl der kanten im graphen schon ) */
             } BRANCH;

typedef struct wurzel {
             unsigned char tiefe;
             unsigned char wurzeltot;
             struct zweig zweige;
	     }  ROOT;

typedef unsigned char vertauschung[2];

typedef signed char QUALITY_1[2*(knoten+1)];
            /* Die Qualitaet eines festen Knotens bei fester Herkunft des
               Baumes */
typedef QUALITY_1 QUALITY_2[knoten+1];
            /* Ein Vektor von Qualitaeten bei fester Herkunft fuer
               verschiedene Knoten */
typedef QUALITY_2 QUALITY[knoten+1];
            /* Ein Vektor aller QUALITY_2er der alle moeglichen Herkuenfte
               abdeckt. Herkunft=0 bedeutet, dass gerade hier mit der
               Nummerierung begonnen wurde und in jede Richtung ein Baum
               geht */

typedef unsigned char QUALITYLAENGE[knoten+1][knoten+1];


/**************************GLOBALE VARIABLEN************************/

unsigned char knotenzahl;
int kantenzahl, maxvalence, codelaenge, border=0, treeborder, regular=0;
int zykelkanten, maxzykelkanten, taille, output, bipartite;
   /* zykelkanten gibt die (momentane) Anzahl der zykelschliessenden Kanten an
      und maxzykelkanten die maximal erlaubte */
int colour[knoten+1]; /* fuer bipartite Graphen */
int anzcolours[2]; /* wieviel gibt es bis jetzt mit Farbe 0 oder 1 */
vertauschung *vert, * combinations[8][8];
int comblaenge[8][8]; /* comblaenge[a][b] gibt die Anzahl der Vertauschungen
                         die man braucht um alle permutationen aufzuzaehlen,
                         die a+b elemente so vertauschen, dass die reihenfolge
                         der ersten a elemente respektiert wird */
int fak[9]; /* Der Maximalgrad wird auch in Normalweiter, baumbau ... auf
               8 beschraenkt */
unsigned int MASK[65];
int nr_of_nodes[knoten+1]; 
                 /* nr_of_nodes[i] ist die Anzahl der Knoten mit Valenz i,*/
                 /* die noch festgelegt werden muessen                    */
                 /* Wird ein Knoten mit valenz i als solcher festgelegt,  */
                 /* d.h. geht die Konstruktion schon beim naechsten Knoten*/
                 /* weiter, so wird nr_of_nodes[i] um 1 verkleinert       */
int orig_non[knoten+1]; /* hier merkt man sich die nr_of_nodes Werte, ohne */
                        /* sie zu aendern */
int stillness[knoten+1]; /* Wieviel neue Knoten mit einer Valenz von >= i */
                         /* werden noch gebraucht ?                       */
                         /* Am Anfang: stillness[1]=Anzahl vertices       */
                         /* Aber wenn alle vertices bereits im Graphen    */
                         /* existieren, ist stillness[1]=0 --- obwohl z.B.*/
                         /* stillness[2]>0 moeglich ist                   */
int recover; /* die Optionen */
unsigned char *globalliste; /* Die Liste der bereits erzeugten Graphen */
int graphenzahl;            /* Die Anzahl der bereits in der Liste enthaltenen
                               Graphen */
unsigned char multiplicity; /* Die Multiplizitaet des Graphen -- festgelegt */
                            /* durch die Multiplizitaet bei Knoten 1        */
unsigned char fertig[knoten+1]; /* Wird am Anfang von minidar belegt. 
                   fertig[i] ist true genau dann, wenn Knoten i fertig 
                   konstruiert ist */
char codename[250], lastgraphname[250];
long int Sicherungszeit;
ROOT roots[knoten+1];
GRAPH lastgraph;  /* Wird fuer recover gebraucht */

int blaetterzahl[knoten+1]; /* die jeweilige anzahl der in den baeumen
                               vorhandenen Blaetter; blaetterzahl[0]
                               enthaelt den wert von g[0][1], wann die
                               baumanalyse durchgefuehrt wurde (default: 0) */
BRANCH **(blaetterkette[knoten+1]);
                            /* Hier werden die Blaetter drin gespeichert, so
                               dass bei der endueberpruefung direkt darauf
                               zugegriffen werden kann. Belegt wird das array
                               bei der analyse durch baumanalyse */

unsigned char * reihenfolge[knoten+1];
            /* rhf[vertex] ist die Liste der zu vertex adjacenten Knoten,
               bei denen die Reihenfolge festgelegt ist. Der Platz hierfuer
               wird am Anfang von main allociert. Das Ende wird durch eine
               0 markiert. Bei maximalvalenz 8 reichen also 9 elemente */

unsigned int drin[knoten+1][2];
                    /* drin[x] ist eine binaere Darstellung der Menge aller
                       Knoten, die zu x adjacent sind und bei denen die
                       Reihenfolge schon festgelegt ist. drin[x] ist also
                       0 genau dann, wenn so verfahren werden soll wie bisher.
                       Val-1 Knoten werden weiterhin extra behandelt, da hier
                       ja nicht nur die Reihenfolge untereinander, sondern
                       sogar die Stellung (am Schluss) festliegt. Leider
                       braucht man 2 long int um ueber 32 Knoten zu
                       behandeln */

int mod=0, rest, splitlevel, splitzaehler=0, noout=1;
int number_of_graphs=0;
int dbz=0;

int debug = 0;

void reset() {
    // set all global variables to what they were at the beginning
    knotenzahl = 0;
    kantenzahl = 0, maxvalence = 0, codelaenge = 0, border=0, treeborder = 0, regular=0;
    zykelkanten = 0, maxzykelkanten = 0, taille = 0, output = 0, bipartite = 0;
    memset(colour, 0, sizeof(colour));
    memset(anzcolours, 0, sizeof(anzcolours));
    memset(comblaenge, 0, sizeof(comblaenge));
    memset(fak, 0, sizeof(fak));
    memset(MASK, 0, sizeof(MASK));
    memset(nr_of_nodes, 0, sizeof(nr_of_nodes));
    memset(orig_non, 0, sizeof(orig_non));
    memset(stillness, 0, sizeof(stillness));
    recover = 0;
    graphenzahl = 0;
    multiplicity = 0;
    memset(fertig, 0, sizeof(fertig));
    memset(codename, 0, sizeof(codename));
    memset(lastgraphname, 0, sizeof(lastgraphname));
    Sicherungszeit = 0;
    memset(roots, 0, sizeof(roots));
    memset(lastgraph, 0, sizeof(lastgraph));
    memset(blaetterzahl, 0, sizeof(blaetterzahl));
    memset(blaetterkette, 0, sizeof(blaetterkette));
    memset(reihenfolge, 0, sizeof(reihenfolge));
    memset(drin, 0, sizeof(drin));
    mod = 0, rest = 0, splitlevel = 0, splitzaehler = 0, noout = 1;
    number_of_graphs = 0;
    dbz = 0;
    vert = NULL;
    memset(combinations, 0, sizeof(combinations));
}

/**************************DIE FUNKTIONEN*******************************/

/*
schreibebaum(BRANCH *zweig, unsigned char biswo)
//BRANCH * zweig;
//unsigned char biswo;
{
int i;
fprintf(stderr,"anfang schreibebaum\n");
schreibezweig(zweig,biswo);
for(i=0; i<zweig->anzahlnachfolger; i++)
     schreibebaum(zweig->successors[i],biswo);
}
*/

/**********************************************************************/

void schreibebild(BILD bild, unsigned char wieweit)

// BILD bild;
// unsigned char wieweit;
{
int i;

for(i=1; i<=wieweit; i++) fprintf(stderr,"%d->%d\n",i,bild[i]);
fflush(0); }

/************************************************************************/
/*
schreibezweig(BRANCH *zweig, unsigned char biswo)

// BRANCH * zweig;
// unsigned char biswo;
{
int i;
fprintf(stderr,"\nAdr.des zweiges: %d \n ",zweig);
fprintf(stderr,"wahrer name des zu bearbeitenden Knotens: %d \n",zweig->wo);
fprintf(stderr,"Blatt: %d \n",zweig->blatt);
fprintf(stderr,"Tiefe: %d \n",zweig->tiefe);
fprintf(stderr,"Zweignummer: %d \n",zweig->zweignummer);
fprintf(stderr,"anzahlnachfolger:%d\n",zweig->anzahlnachfolger);
fprintf(stderr,"anzahltot %d  \n ",zweig->anzahltot );
fprintf(stderr,"Adr.des nachfolgerarrays: %d \n ",zweig->successors);
for (i=0;i<zweig->anzahlnachfolger; i++)
    fprintf(stderr,"Adr.von nachf %d: %d  \n ",i,zweig->successors[i] );
for (i=0;i<zweig->anzahlnachfolger; i++)
    fprintf(stderr,"tot[%d]: %d  \n ",i,zweig->tot[i] );
fprintf(stderr,"Adr.von vorg: %d \n ",zweig->vorgaenger );
fprintf(stderr,"Rootnummer: %d \n",zweig->rootnumber);
schreibebild(zweig->bild,biswo);
fprintf(stderr,"\n\n");
fflush(0);
}
*/

/**************************************************************************/

void schreibegraph(GRAPH g)

//GRAPH g;
{
 int x,y;
fprintf(stderr,"\n\n ");

fprintf(stderr,"|%2d",g[0][0]);

for(x=1; (x <= g[0][0])&&(x<=24); x++)  fprintf(stderr,"|%2d",x); fprintf(stderr,"|\n");

fprintf(stderr," ");

for(x=0; (x <= g[0][0])&&(x<=24); x++) fprintf(stderr,"|==");    fprintf(stderr,"|\n");

for(x=0; x < maxvalence; x++)
  {
   fprintf(stderr," |  ");
   for(y=1; (y<=g[0][0])&&(y<=24); y++)
       if (g[y][x] ==leer) fprintf(stderr,"|  "); else fprintf(stderr,"|%2d",g[y][x]);
       fprintf(stderr,"|\n");
  }

if (g[0][0]>24) 
{
fprintf(stderr,"\n");

fprintf(stderr,"    ");

for(x=25; (x <= g[0][0])&&(x<=48); x++)  fprintf(stderr,"|%2d",x); fprintf(stderr,"|\n");

fprintf(stderr,"    ");

for(x=25; (x <= g[0][0])&&(x<=48); x++) fprintf(stderr,"|==");    fprintf(stderr,"|\n");

for(x=0; x < maxvalence; x++)
  {
   fprintf(stderr,"    ");
   for(y=25; (y <= g[0][0])&&(y<=48); y++)
       if (g[y][x] ==leer) fprintf(stderr,"|  "); else fprintf(stderr,"|%2d",g[y][x]);
       fprintf(stderr,"|\n");
  }
}

if (g[0][0]>48) 
{
fprintf(stderr,"\n");

fprintf(stderr,"    ");

for(x=49; x<=g[0][0]; x++)  fprintf(stderr,"|%2d",x); fprintf(stderr,"|\n");

fprintf(stderr,"    ");

for(x=49; x <= g[0][0]; x++) fprintf(stderr,"|==");    fprintf(stderr,"|\n");

for(x=0; x < maxvalence; x++)
  {
   fprintf(stderr,"    ");
   for(y=49; y<=g[0][0]; y++)
       if (g[y][x] ==leer) fprintf(stderr,"|  "); else fprintf(stderr,"|%2d",g[y][x]);
       fprintf(stderr,"|\n");
  }
}
fflush(0);}


/************************KONSTRMENGE**************************************/


void konstrmenge(GRAPH graph, unsigned char vertex, unsigned char max,
                 unsigned char tiefe, unsigned int *menge, unsigned int *vmg)

/* hier werden knoten mit abstand tiefe der Menge hinzugefuegt */

//GRAPH graph;
//unsigned char vertex, max, tiefe;
//unsigned int *menge;
//unsigned int *vmg; /* vergleichsmenge */

{  
int i;

for (i=0; graph[vertex][i]!=leer; i++)
  { if (graph[vertex][i]<=32)
	  { menge[0] = menge[0] | MASK[graph[vertex][i]];
         if ( (tiefe < max) && !(menge[0] & vmg[0]) ) 
                konstrmenge(graph,graph[vertex][i],max,tiefe+1,menge,vmg);
	  }
        else
	  { menge[1] = menge[1] | MASK[graph[vertex][i]];
         if ( (tiefe < max) && !(menge[1] & vmg[1]) ) 
                konstrmenge(graph,graph[vertex][i],max,tiefe+1,menge,vmg);
	  }
  }
} 

/*******************************ABSTANDOK******************************/

int abstandok(GRAPH graph, int v1, int v2, int taille)

//GRAPH graph;
//int v1, v2;
//int taille;

/* ueberprueft, ob der Abstand von v1 und v2 in graph <= taille-1 ist */
/* das geschieht, indem die Menge m1 aller Knoten, die einen Abstand */
/* von (taille-2)/2 von v1 haben  ( und analog m2 bezueglich v2 ) konstruiert*/
/* werden und der schnitt daraufhin ueberprueft wird, ob er leer ist */

{

unsigned int m1[2], m2[2];

if (taille<=3) return(1);
/* Dass fuer girth==3 der returnwert immer ok ist, wird dadurch sichergestellt,
   dass die maximale multiplizitaet immer bei 1->2 auftritt -- in konstrukt
   werden also deshalb keine doppelkanten eingefuegt, weil das in v1_konstrukt
   nicht gemacht wurde */



if (v1<=32) { m1[0]=MASK[v1]; m1[1]=0; }
       else { m1[1]=MASK[v1]; m1[0]=0; }



if (v2<=32) { m2[0]=MASK[v2]; m2[1]=0; }
       else { m2[1]=MASK[v2]; m2[0]=0; }

konstrmenge(graph,v1,((taille-1)/2),1,m1,m2);

if (!((m1[0] & m2[0]) || (m1[1] & m2[1]) ) )
konstrmenge(graph,v2,((taille-2)/2),1,m2,m1);

return( !((m1[0] & m2[0]) || (m1[1] & m2[1]) ) );

}


/**********************************CODIERE********************************/

void codiere(GRAPH graph, unsigned char *code, ADJAZENZ adj)

//GRAPH graph;
//unsigned char *code;
//ADJAZENZ adj;

{
int i,j,codestelle;

/* Codierung: Erstes byte: knotenzahl danach die zu knoten 1 adjazenten
   ALLE GROESSEREN vertices, beendet durch 0, danach die zu knoten 2 
   adjazenten GROESSEREN .... jeweils beendet durch 0. Da der letzte Knoten
   nie zu groesseren adjazent sein kann, braucht seine (immer leere) Liste 
   nicht durch 0 abgeschlossen werden*/


code[0]=knotenzahl;

for (codestelle=1; codestelle<=adj[1]; codestelle++) 
                         code[codestelle]=graph[1][codestelle-1]; 
/* 0 ist das Trennzeichen fuer den naechsten Code */
/* Gesamtcodelaenge: Kantenzahl+knotenzahl      */

for (j=2; j<=graph[0][0]; j++)
{ code[codestelle]=0; codestelle++;
  for (i=1; i<adj[j]; i++) 
                     if (graph[j][i]>j) 
                            { code[codestelle]=graph[j][i]; codestelle++; }
}

}

/********************************WEGSPEICHERN******************************/

void wegspeichern(unsigned char *liste, int graphenzahl,
                  GRAPH graph, unsigned char vertex)

// unsigned char *liste;
// int graphenzahl;
// GRAPH graph;
// unsigned char vertex;
{
FILE *fil;

/* Zuerst muss sich noch der vertex gemerkt werden, bei dem man war */
graph[0][knoten]=vertex;

/*fil=fopen(lastgraphname,"w");
fwrite(graph,sizeof(GRAPH),1,fil);
fclose(fil);*/

graph[0][knoten]=0;

if (output)
  { fil=stdout;
    // fwrite(liste,codelaenge,graphenzahl,fil);
  }
 else
   if (!noout)
  {
fil=fopen(codename,"a");
fwrite(liste,codelaenge,graphenzahl,fil);
fclose(fil);
}
}


/**********************************AUFSCHREIBEN*****************************/

void aufschreiben(GRAPH g, ADJAZENZ adj)

//GRAPH g;
//ADJAZENZ adj;

/*** Nimmt den Graphen und schreibt ihn vorerst in die globale Liste */
/* wenn die voll ist, wird sie weggesichert. */

{

number_of_graphs++;
if (noout) return;
codiere(g,globalliste+(graphenzahl*codelaenge),adj);

graphenzahl++;
if (graphenzahl==listenlaenge) 
    {     wegspeichern(globalliste,graphenzahl,g,knotenzahl+1); graphenzahl=0;

	     }


}

/******************baumeinfuegen***********************************/

void baumeinfuegen(QUALITY_2 qual, unsigned char *liste,
                   int neuer, unsigned char *laenge)

//QUALITY_2 qual; /* Die fuer diese herkunft geltende Liste der Qualitaeten */
//unsigned char * liste; /* die liste in die der neue nachfolger einsortiert 
//                          werden muss */
//int neuer;   /* der neue nachfolger */
//unsigned char * laenge; /* die laengen der qualitaetslisten */

{
signed char *q1, *q2;
int j,k,co,grenze;

/*fprintf(stderr,"anf be neuer: %d\n",neuer);
for(j=0; j<8;j++) fprintf(stderr," %d ",liste[j]); fprintf(stderr,"be davor\n");*/

for (j=0;liste[j];j++); /* sucht Ende der Liste */
liste[j+1]=0;
co= 1;
q1=qual[neuer];


for ( j-- ;(co > 0) && (j>=0); j-- )
   {
     if (laenge[liste[j]]>laenge[neuer]) grenze=laenge[neuer];
                             else grenze=laenge[liste[j]];
        q2=qual[liste[j]];
        for (k=0; (q1[k]==q2[k]) && (k<=grenze); k++); 
                   /* sucht Unterschied oder Ende "<=" weil das endsymbol
                      in der laenge nicht mitgezaehlt wird */
        co = q1[k]-q2[k];
        if (co>0) liste[j+1]=liste[j]; else liste[j+1]=neuer;
       }
  if ((j<0) && (co>0)) liste[0]=neuer;

/*for(j=0; j<8;j++) fprintf(stderr," %d ",liste[j]); fprintf(stderr,"be danach\n");*/

}


/*******************qualitaetbelegen**********************************/

void qualitaetbelegen(unsigned char vertex, unsigned char valence, QUALITY quality,
                      signed char *zubelegen, int anzahl, unsigned char *reihenfolge, 
                      unsigned char *hlaenge, signed char *neuelaenge,
                      unsigned char luecke, unsigned char anzahlnull )

//unsigned char vertex;
//unsigned char valence; /* adj[vertex] */
//QUALITY quality;
//signed char * zubelegen;
//int anzahl;
//unsigned char * reihenfolge;
//unsigned char * hlaenge; /* das laengensegment, das dieser herkunft 
//                           entspricht */
//signed char * neuelaenge; /* das neu zu belegende laengensegment */
//unsigned char luecke; /* welche qualitaet soll nicht in die neue hineingemerged
//                         werden, da sie der herkunft entspricht */
//unsigned char anzahlnull; /* wieviel nuller muessen in das erste niveau
//                             eingefuegt werden ? */

{
int wieweit[knoten+1]; /* wieweit habe ich bei dem i-ten Element die
                          Qualitaet schon abgeschrieben ? */
int i, zaehler, ausfall, biswo;
unsigned char which, ll;
signed char *qual;

/*fprintf(stderr,"anfang qb. vertex: %d luecke: %d anzahl: %d\n",vertex,luecke,anzahl);
for(i=0; i<8;i++) fprintf(stderr," %d ",reihenfolge[i]); fprintf(stderr,"\n");*/

if (luecke==leer) biswo=anzahl; else biswo=anzahl+1;

for (i=0; i<biswo; i++) wieweit[i]=0;

zubelegen[0]=valence-1; zubelegen[1]= -1;
zaehler=2;

ausfall=0;
/* Im ersten Niveau muessen auch die Valenz 1 Knoten notiert werden, deshalb
   eine Sonderbehandlung: */
for (i=0; i<biswo; i++)
          if (reihenfolge[i] != luecke)
           {
           which=reihenfolge[i];
           ll=hlaenge[which];
           if (wieweit[i]<ll)
             { qual=quality[vertex][which];
               for ( ;qual[wieweit[i]]!= -1; (wieweit[i])++)
                   { zubelegen[zaehler]=qual[wieweit[i]];
                     zaehler++; 
		   }
               if (wieweit[i]==ll) ausfall++;
               wieweit[i]++;
             }
            } /* Die hier zu ende gehende schleife schreibt ein ganzes niveau
                 des baumes */
            for(i=0;i<anzahlnull;i++) { zubelegen[zaehler]=0; zaehler++; }
            zubelegen[zaehler]= -1;
            zaehler++;


for(; ausfall<anzahl; )
       {
        for (i=0; i<biswo; i++)
          if (reihenfolge[i] != luecke)
           {
           which=reihenfolge[i];
           ll=hlaenge[which];
           if (wieweit[i]<ll)
             { qual=quality[vertex][which];
               for ( ;qual[wieweit[i]]!= -1; (wieweit[i])++)
                   { zubelegen[zaehler]=qual[wieweit[i]];
                     zaehler++; }
               if (wieweit[i]==ll) ausfall++;
               wieweit[i]++;
             }
            } /* Die hier zu ende gehende schleife schreibt ein ganzes niveau
                 des baumes */
            zubelegen[zaehler]= -1;
            zaehler++;
       }
zubelegen[zaehler]= -2;
(*neuelaenge)=zaehler-1;
/*fprintf(stderr,"ende qb\n");*/
}

/******************bearbeite**************************************/

void bearbeite(GRAPH graph, ADJAZENZ adj, unsigned char *marks, unsigned char vertex,
               unsigned char *liste, unsigned char *anzahl, unsigned char neuetiefe,
               unsigned char *wieviel, QUALITY quality, QUALITYLAENGE laenge)

//GRAPH graph;
//ADJAZENZ adj;
//unsigned char * marks;
//unsigned char vertex; /* welchen vertex bearbeite ich */
//unsigned char * liste;/* die liste der in der naechsten tiefe zu bearbeitenden
//                         Knoten */
//unsigned char * anzahl; /* die anzahl der in der naechsten tiefe zu 
//                           bearbeitenden Knoten */
//unsigned char neuetiefe;
//unsigned char * wieviel; /* wieviel baeume haengen schon an diesem Knoten */
//QUALITY quality;
//QUALITYLAENGE laenge;

{

static unsigned char localcount, herkunft=0, i, ziel, anzahlnull;

/* anzahlnull ist die anzahl der Nullqualitaeten, die eingefuegt werden 
   muessen */

/*fprintf(stderr,"anfang bearb vertex: %d wieviel: %d\n",vertex,*wieviel);*/

localcount=0;

for (i=0;i<adj[vertex];i++)
         { ziel=graph[vertex][i];
           if (laenge[vertex][ziel])
                { localcount++;
                  if (((ziel<=32) && !(MASK[ziel] & drin[vertex][0])) ||
                      ((ziel>32) && !(MASK[ziel] & drin[vertex][1])))
                           /* d.h.: Da haengt ein neuer Baum dran */
                   { 
                   if (ziel<=32) drin[vertex][0] |= MASK[ziel];
                          else drin[vertex][1] |= MASK[ziel];
                   (*wieviel)++;
                   baumeinfuegen(quality[vertex],reihenfolge[vertex],
                                                     ziel,laenge[vertex]);
                   }
                }
            else if (adj[ziel]!=1) herkunft=ziel; 

          }


/* Jetzt weiss ich genau, wieviel baumanhaenger ich habe: *wieviel --- und
   wieviel in meiner Reihenfolgetabelle stehen: localcount.
   Also habe ich ((*wieviel) - localcount) adjacente valenz 1 Knoten */


anzahlnull=(*wieviel)-localcount;

if ((*wieviel)==adj[vertex]-1) /* der baum geht eindeutig weiter */
            {
              if (!(marks[herkunft]==neuetiefe))
                 { 
                   liste[(*anzahl)]=herkunft;
                   (*anzahl)++;
                   marks[herkunft]=neuetiefe;
                 }
                 qualitaetbelegen(vertex,adj[vertex],quality,
                      quality[herkunft][vertex],localcount,reihenfolge[vertex],
                      laenge[vertex],&(laenge[herkunft][vertex]),
                      leer, anzahlnull );
             }

else
if ((*wieviel)==adj[vertex]) /* Es ist durch und durch ein Baum --
                                man kommt von allen Seiten */
     { 
     localcount--; /* Einer aus der Liste wird ja jeweils nicht mit in
                                 die Qualitaet einbezogen */
       for (i=0; i<adj[vertex]; i++)
         { herkunft=graph[vertex][i];
           if ((!laenge[herkunft][vertex]) && (adj[herkunft]>1))
            { 
              if (!(marks[herkunft]==neuetiefe))
                 { liste[*anzahl]=herkunft;
                   (*anzahl)++;
                   marks[herkunft]=neuetiefe;
                 }
                 qualitaetbelegen(vertex,adj[vertex],quality,
                     quality[herkunft][vertex],localcount,reihenfolge[vertex],
                     laenge[vertex],&(laenge[herkunft][vertex]),
                     herkunft, anzahlnull);

             }
           }
      }

}
           
/***********************preprocessing*******************************/

void preprocessing(GRAPH graph, ADJAZENZ adj)

//GRAPH graph;
//ADJAZENZ adj;

{
int i,j;
unsigned char herkunft, tiefe;
unsigned char  marks[knoten+1], wieviel[knoten+1], anzahl[knoten+1];
/* marks[i]=tiefe wenn Knoten i in der Tiefe tiefe bearbeitet werden muss.
   wieviel[i] gibt an, wieviel baumartige fortsaetze knoten i schon hat.
   anzahl[i] gibt an, wieviel Knoten es bis jetzt in Tiefe i gibt.
   Diese Variablen dienen der Organisation der "Breitensuche" von den
   Valenz-1-Knoten aus */
unsigned char liste[knoten+1][knoten+1];/*Die eigentliche Breitensucheliste*/
QUALITY quality; /* das 3-dim array in dem die Baumqualitaeten gespeichert
                    werden */
QUALITYLAENGE laenge;    /* wie lang der Qualitaets-
                            identifier von herkunft nach vertex schon ist,
                            steht in laenge[herkunft][vertex] */


for (i=1; i<=graph[0][0]; i++) { marks[i]=anzahl[i]=drin[i][0]=drin[i][1]=0;
                                 wieviel[i]=reihenfolge[i][0]=0;
                                 for(j=0;j<adj[i]; j++) 
                                    laenge[i][graph[i][j]]=0; }


for (i=1; i<=graph[0][0]; i++)
    if (adj[i]==1)
       { herkunft=graph[i][0]; 
         wieviel[herkunft]++;
               /* Sonderbehandlung fuer Valenz 1 Knoten -- nicht als 
                  dranhaengende Baeume bearbeiten -- nur "wieviel" 
                  hochzaehlen, da in "bearbeite" solche Knoten nicht
                  beachtet werden */
         if (!marks[herkunft]) /* d.h. die herkunft wird zum ersten mal
                                  erreicht */
             { liste[1][anzahl[1]]=herkunft;
               marks[herkunft]=1; 
               anzahl[1]++;
             } /* Neu in die Liste der zu bearbeitenden Knoten aufnehmen */
        }



for ( tiefe=1; anzahl[tiefe]!=0; tiefe++ )
  { 
    for (i=0; i<anzahl[tiefe]; i++)
    bearbeite(graph,adj,marks,liste[tiefe][i],liste[tiefe+1],anzahl+tiefe+1,
                     tiefe+1,&(wieviel[liste[tiefe][i]]),quality,laenge);
   }

for (i=1; i<=graph[0][0]; i++) if (reihenfolge[i][1]==0) 
                                     drin[i][0]=drin[i][1]=0; 
    /* Wenn nur ein Baum dranhaengt ist es besser diese Verbindung als normale
       einser-Adjazenz zu behandeln */

}


/****************************ENDSUCHE**********************************/

void endsuche(BRANCH *zweig, BRANCH **blaetterkette, int *zaehler)

//BRANCH *zweig;
//BRANCH **blaetterkette;
//int *zaehler; /* das wievielte element der kette muss beschrieben werden */

{ int i;


if ( zweig->blatt ) 
    { blaetterkette[*zaehler]=zweig; (*zaehler)++; }

else 
   for (i=0; i<zweig->anzahlnachfolger; i++)
     if ( (zweig->successors[i] != nil) && (!zweig->tot[i]) )
                          endsuche(zweig->successors[i],blaetterkette,zaehler);
}


/**************************BAUMANALYSE**************************************/

void baumanalyse(unsigned char bis_wo)

//unsigned char bis_wo; /* gibt an, bis zu welcher zahl eine analyse noetig
//                         ist. Im allgemeinen wird graph[0][0] uebergeben */

{
int i,zaehler;

for (i=1; i<=bis_wo; i++)
    { 
      if ((roots[i].tiefe) && (!(roots[i].wurzeltot)))
	  { zaehler=0;
          blaetterkette[i]=(BRANCH **)malloc(sizeof(BRANCH *)*blaetterzahl[i]);
                     endsuche(&(roots[i].zweige),blaetterkette[i],&zaehler);
	  }
     }
}

/**************************ZWEIGENTFERNEN*****************************/

void zweigentfernen(BRANCH *zweig)

//BRANCH *zweig;

{ int i;

  if (zweig->blatt)  blaetterzahl[zweig->rootnumber]--;
    else
     {  for (i=0; i<zweig->anzahlnachfolger; i++)
          if (zweig->successors[i] != nil) 
                        zweigentfernen(zweig->successors[i]);
        free(zweig->tot);
        free(zweig->successors); }
  free(zweig); 
}
 


/****************************ROOTENTFERNEN********************************/

void rootentfernen(unsigned char i)

//unsigned char i;

/* entfernt die i-te wurzel (genauer: entfernt alles, was daran haengt)
   und initialisiert die werte der wurzel */

{
int j;
BRANCH *zweig;

zweig= &(roots[i].zweige);
for (j=0; j<zweig->anzahlnachfolger; j++) 
      if (zweig->successors[j] != nil)
           zweigentfernen( zweig->successors[j] );
free(zweig->successors);
free(zweig->tot);

roots[i].wurzeltot=roots[i].tiefe=zweig->anzahlnachfolger=0;
zweig->anzahltot=blaetterzahl[i]=0;
}



/**********************ZWEIGAUFRAEUMEN**********************************/


BRANCH *zweigaufraeumen(BRANCH *zweig, unsigned char k)

//BRANCH *zweig;
//unsigned char k;

/* ich komme hier nur hin, wenn diese abzweigung nicht schon in einem
   schritt vor k totgesetzt wurde */

{ 
int j;

if (zweig->tiefe == k) { 
                       if (zweig->blatt) blaetterzahl[zweig->rootnumber]--;
                       else
		       {
                         for (j=0;j<zweig->anzahlnachfolger;j++)
                                  if (zweig->successors[j] != nil) 
                                       zweigentfernen(zweig->successors[j]);
                         free(zweig->successors);
                         free(zweig->tot);
		       }
                       free(zweig);
                       return(nil); }
else /* d.h. es ist ein aelterer zweig */
{
if (!zweig->blatt)
   {
   for (j=0;j<zweig->anzahlnachfolger;j++)
      {
      if (zweig->tot[j] == k) { zweig->tot[j]=0; (zweig->anzahltot)--; }
      if ((zweig->successors[j] != nil) && (zweig->tot[j]==0) )
           zweig->successors[j]=zweigaufraeumen(zweig->successors[j],k);
      }
    if (zweig->anzahltot == 0) zweig->blatt = 1;
               /* Moeglicherweise ein neues Blatt */
    for (j=0;(j<zweig->anzahlnachfolger) && zweig->blatt;j++)
                if (zweig->successors[j]!=nil) zweig->blatt=0;
    if (zweig->blatt) {
                        blaetterzahl[zweig->rootnumber]++;
                        zweig->anzahlnachfolger=0;
                        free(zweig->successors);
                        zweig->successors=nil;
                        free(zweig->tot);
                        zweig->tot=nil; }
   }
return(zweig);
}/* ende else */
}




/****************************BAUMAUFRAEUMEN******************************/


void baumaufraeumen(unsigned char k, unsigned char vertexzahl)

//unsigned char k;
//unsigned char vertexzahl;

/* entfernt alle im schritt k gesetzten Blaetter und totmarken 
   ueberprueft dabei alle an den roots 1..vertexzahl haengenden Baeume .
   blaetter und totmarken mit einem index > k existieren zu diesem Zeitpunkt
   schon nicht mehr im baum */

{

int i,j;


for (i=1; i<=vertexzahl; i++)
    { if (roots[i].tiefe)
         {
            if (roots[i].tiefe == k) rootentfernen(i);
              else /* d.h.: es ist eine schon aeltere wurzel */
                  { 
                     if (roots[i].wurzeltot==k) roots[i].wurzeltot=0;
                     if (!(roots[i].wurzeltot))
		       {
                       for (j=0; j<roots[i].zweige.anzahlnachfolger; j++)
                          { if (roots[i].zweige.tot[j]==k) 
                              { roots[i].zweige.tot[j]=0;
                                (roots[i].zweige.anzahltot)--;
			      }
                            if (!roots[i].zweige.tot[j])
                                 roots[i].zweige.successors[j]=
                              zweigaufraeumen(roots[i].zweige.successors[j],k);
	 		   
/* wenn dieser ast schon vorher tot war, wurde er auch im k-ten schritt
   nicht mehr benutzt */
			  }
		     } /* ende not wurzeltot */
                  } /* ende aeltere wurzel */
          } /* ende roots[i].tiefe -- d.h. schon initialisiert */
      }

} /* ende baumaufraeumen */

/*********************************ordnen******************************/

void ordnen(unsigned char *tr, int anzahl)
/* ordnet einfach nur die naechsten anzahl char nach dem Zeiger tr. */
/* Bubblesort */

//unsigned char* tr;
//int anzahl;

{

unsigned char puffer, getauscht;
int i,j;

getauscht=1;
for (i=anzahl-1; (i>0) && getauscht ; i--)
  { getauscht=0;
    for (j=0; j<i; j++) if (tr[j+1]<tr[j]) { getauscht=1;
                                             puffer=tr[j];
                                             tr[j]=tr[j+1];
                                             tr[j+1]=puffer;
					   }
  }


}

/****************************COMPARE************************************/

int compare (unsigned char bild[knoten+1], unsigned char *tr, unsigned char *br)

//unsigned char bild[knoten+1];  /* die Bilder der Knoten */
//unsigned char *tr;             /* die Reihe in der die Bilder eingesetzt */
//                               /* werden muessen                         */
//unsigned char *br;             /* die Reihe mit der die ersetzte verglichen */
//                               /* wird */

/* vergleicht die Reihe tauschreihe, in der alle Elemente durch ihre Bilder */
/* ersetzt werden mit bestreihe. Gibt etwas negatives zur"uck wenn tr<br, 0 */
/* wenn tr=br und etwas positives sonst */

{

unsigned char test[knoten+1];
int i;

test[0]=bild[tr[0]];
for (i=1; tr[i]!=leer; i++) test[i]=bild[tr[i]];
test[i]=leer;
ordnen(test,i); /* i ist als index eins zu gross, als Anzahl aber passend */
for (i=0; (test[i]==br[i]) && br[i]!=leer; i++);

return((int)(test[i]-br[i]));

}



/***********************NORMALWEITER********************************/


void normalweiter(GRAPH graph, ADJAZENZ adj, BILD bild, unsigned char wo,
                  unsigned char naechstenummer, int *abbruch)

//GRAPH graph;
//ADJAZENZ adj;
//BILD bild;
//unsigned char wo; /* wahrer namen des hier zu behandelnden knotens */
//unsigned char naechstenummer; /* naechste zu vergebende nummer */
//int *abbruch;

{
int i,j,k,suche,co;
unsigned char merke;
int vertex, zw1, zw2, puffer;
unsigned char *(candidates[3]); /* Die Loecherkandidaten, fuer die kleinste 
          nummer werden in dem vector, auf den candidates[0],... zeigt, 
          gespeichert und entsprechend weiter. Bei Auslegung auf Maximalvalenz
          8 reichen 3 Zeiger -- ueberlegen*/
unsigned char kandidatenspeicher[9][9]; /* In Kandidatenspeicher[i][0] wird 
          die Anzahl der Kandidaten mit multiplizitaet[i] gespeichert und in 
          den folgenden elementen die kandidaten selbst */
 unsigned char merkebild[9];
unsigned char einserspeicher[9]; /* Hier werden die knoten mit grad 1 
          gespeichert -- sie bekommen immer die letzten Zahlen und die werden
          auch nicht vertauscht */
unsigned char rhfspeicher[9]; /* Hier werden die knoten gespeichert, deren
          Reihenfolge durch die Baumstruktur festgelegt ist */

int maxmultiplicity, numberofcand;
          /* maximale multiplizitaet des bearbeiteten Knotens und die Anzahl
             der verschiedenen multiplizitaeten */
unsigned char ziel, a, fixed; /* fixed gibt an, wieviel der Kanten schon in der
                              Reihenfolge festgelegt sind */
unsigned int hierdrin[2];
vertauschung * transpo;


if (adj[wo]>0)
{
   hierdrin[0]=hierdrin[1]=0;
   maxmultiplicity=1;
   for (i=0; i<=maxvalence; i++) kandidatenspeicher[i][0]=0;
   einserspeicher[0]=rhfspeicher[0]=0;
   vertex=bild[wo];
   suche= -1;

if (drin[wo][0] | drin[wo][1]) 
                  /*der else teil ist fast identisch -- nur die 
                    drin-spezifischen Abfragen spart man sich halt */
{
   for (i=0; i<adj[wo]; )
     { ziel=graph[wo][i];
       if (bild[ziel]==unbelegt)
        { 
          if (adj[ziel]==1)
             { (einserspeicher[0])++;
               einserspeicher[einserspeicher[0]]=ziel;
               i++;
               suche=1;}
         else
           if ( ((ziel<=32) && (MASK[ziel] & drin[wo][0])) ||
                ((ziel>32) && (MASK[ziel] & drin[wo][1])) )
             { 
               if (ziel<=32) hierdrin[0] |= MASK[ziel];
                  else hierdrin[1] |= MASK[ziel];
               i++;
               suche=1;}
         else
          { 
            suche=i; merke=graph[wo][i]; i++;
            for( ;(graph[wo][i]==merke) && (i<adj[wo]); i++); 
            suche=i-suche; kandidatenspeicher[suche][0]++; 
            kandidatenspeicher[suche][kandidatenspeicher[suche][0]]=merke; 
            if (suche>maxmultiplicity) maxmultiplicity=suche; }
	}
         else i++; }


if (suche >=0) /* d.h. es gibt mindestens ein loch */
  {
     numberofcand=0;
      candidates[0]=candidates[1]=candidates[2]=kandidatenspeicher[0]; 
                         /* Diese werte werden nicht notwendig belegt, 
                         muessen aber fuer die Abfrage in der for-schleife 
                         den Defaultwert 0 beim Element c[i][0] haben */
fixed=0;
 for (j=0;(ziel=reihenfolge[wo][j]);j++)
    if( ((ziel<=32) && (MASK[ziel] & hierdrin[0])) ||
        ((ziel>32) && (MASK[ziel] & hierdrin[1])) )
      { (kandidatenspeicher[1][0])++; fixed++;
        kandidatenspeicher[1][kandidatenspeicher[1][0]]=ziel; }

    for (i=maxmultiplicity; i>0; i--)
      { if (kandidatenspeicher[i][0]) 
                    { candidates[numberofcand]=kandidatenspeicher[i];
                      numberofcand++;
                      for (j=1; j<=kandidatenspeicher[i][0]; j++) 
                              {
                              bild[kandidatenspeicher[i][j]]=naechstenummer;
                              naechstenummer++; }
		    }
      }

for (j=1; j<=einserspeicher[0]; j++)
    { bild[einserspeicher[j]]=naechstenummer;
      naechstenummer++;
    }

if (numberofcand<3) { candidates[2]=candidates[numberofcand-1];
                      candidates[numberofcand-1]=kandidatenspeicher[0]; }

/* das dient dazu, die Val1-knoten und die mit festgelegter reihenfolge
   immer an der gleichen Stelle zu haben -- in der 2-er Schleife */

    co=compare(bild,graph[wo],graph[bild[wo]]);


    if (co<0) *abbruch=1;
    else 
 if ((vertex<graph[0][0])&&(co==0))
    {
    for (suche=1; bild[suche]!=vertex+1; suche++); 
                                 /*sucht den Knoten mit bild vertex+1 */
    a=candidates[2][0]-fixed; transpo=combinations[a][fixed];

    for (j=0; (j<fak[candidates[1][0]]) && !(*abbruch); j++)
    {
    if (j)
    { zw1=candidates[1][vert[j][0]];
    zw2=candidates[1][vert[j][1]];
    puffer=bild[zw1];
    bild[zw1] = bild[zw2];
    bild[zw2]=puffer;
    if (suche==zw1) suche=zw2; else { if (suche==zw2) suche=zw1; } 
                 /* hat sich der Wert des naechsten Knotens veraendert ? */
    }

    for (k=0; (k<fak[candidates[0][0]]) && !(*abbruch); k++)
    { 
    if (k)
    { zw1=candidates[0][vert[k][0]];
    zw2=candidates[0][vert[k][1]];
    puffer=bild[zw1];
    bild[zw1] = bild[zw2];
    bild[zw2]=puffer;
    if (suche==zw1) suche=zw2; else { if (suche==zw2) suche=zw1; } 
                 /* hat sich der Wert des naechsten Knotens veraendert ? */
    }

    for (i=0; (i<comblaenge[a][fixed]) && !(*abbruch); i++)
    {
    if (i)
    { zw1=candidates[2][transpo[i][0]];
      zw2=candidates[2][transpo[i][1]];
      candidates[2][transpo[i][0]]=zw2;
      candidates[2][transpo[i][1]]=zw1;
      for (puffer=1;puffer<=candidates[2][0];puffer++)
	bild[candidates[2][puffer]]=merkebild[puffer];
      if (suche==zw1) suche=zw2; else { if (suche==zw2) suche=zw1; } 
                /* hat sich der Wert des naechsten Knotens veraendert ? */
    }
    else {
      for (puffer=1;puffer<=candidates[2][0];puffer++)
	{ merkebild[puffer]=bild[candidates[2][puffer]]; }
    }

    normalweiter(graph,adj,bild,suche,naechstenummer,abbruch); 
    } /* Ende for [0][0]*/
    } /* Ende for [1][0]*/
    } /* Ende for [2][0]*/

   } /* Ende if vertex<.....*/

    for (j=1; j<=candidates[2][0]; j++) bild[candidates[2][j]]=unbelegt;
    for (j=1; j<=candidates[1][0]; j++) bild[candidates[1][j]]=unbelegt;
    for (j=1; j<=candidates[0][0]; j++) bild[candidates[0][j]]=unbelegt;
    for (j=1; j<=einserspeicher[0]; j++) bild[einserspeicher[j]]=unbelegt;



   }/* Ende es gibt ein Loch */
   else /* d.h. kein loch */
     { 
     co=compare(bild,graph[wo],graph[bild[wo]]);
     if (co<0) *abbruch=1;
     else
     if ((vertex<graph[0][0])&&(co==0))
       {
        for (suche=1; bild[suche]!=vertex+1; suche++); 
                               /*sucht den Knoten mit bild vertex+1 */
        normalweiter(graph,adj,bild,suche,naechstenummer,abbruch); }
     }
}





else /* d.h. not drin[vertex] */
{  for (i=0; i<adj[wo]; )
     { if (bild[graph[wo][i]]==unbelegt)
        { if (adj[graph[wo][i]]==1)
             { (einserspeicher[0])++;
               einserspeicher[einserspeicher[0]]=graph[wo][i];
               i++;
               suche=1;}
         else
          { suche=i; merke=graph[wo][i]; i++;
            for( ;(graph[wo][i]==merke) && (i<adj[wo]); i++); 
            suche=i-suche; kandidatenspeicher[suche][0]++; 
            kandidatenspeicher[suche][kandidatenspeicher[suche][0]]=merke; 
            if (suche>maxmultiplicity) maxmultiplicity=suche; }
	}
         else i++; }

if (suche >=0) /* d.h. es gibt mindestens ein loch */
  {
     numberofcand=0;
      candidates[0]=candidates[1]=candidates[2]=kandidatenspeicher[0]; 
                         /* Diese werte werden nicht notwendig belegt, 
                         muessen aber fuer die Abfrage in der for-schleife 
                         den Defaultwert 0 beim Element c[i][0] haben */
    for (i=maxmultiplicity; i>0; i--)
      { if (kandidatenspeicher[i][0]) 
                    { candidates[numberofcand]=kandidatenspeicher[i];
                      numberofcand++;
                      for (j=1; j<=kandidatenspeicher[i][0]; j++) 
                              {
                              bild[kandidatenspeicher[i][j]]=naechstenummer;
                              naechstenummer++; }
		    }
      }

for (j=1; j<=einserspeicher[0]; j++)
    { bild[einserspeicher[j]]=naechstenummer;
      naechstenummer++;
    }

    co=compare(bild,graph[wo],graph[bild[wo]]);

    if (co<0) *abbruch=1;
    else 
 if ((vertex<graph[0][0])&&(co==0))
    {
    for (suche=1; bild[suche]!=vertex+1; suche++); 
                                 /*sucht den Knoten mit bild vertex+1 */
    for (i=0; (i<fak[candidates[2][0]]) && !(*abbruch); i++)
    {
    if (i)
    { zw1=candidates[2][vert[i][0]];
    zw2=candidates[2][vert[i][1]];
    puffer=bild[zw1];
    bild[zw1] = bild[zw2];
    bild[zw2]=puffer;
    if (suche==zw1) suche=zw2; else { if (suche==zw2) suche=zw1; } 
                /* hat sich der Wert des naechsten Knotens veraendert ? */
    }

    for (j=0; (j<fak[candidates[1][0]]) && !(*abbruch); j++)
    {
    if (j)
    { zw1=candidates[1][vert[j][0]];
    zw2=candidates[1][vert[j][1]];
    puffer=bild[zw1];
    bild[zw1] = bild[zw2];
    bild[zw2]=puffer;
    if (suche==zw1) suche=zw2; else { if (suche==zw2) suche=zw1; } 
                 /* hat sich der Wert des naechsten Knotens veraendert ? */
    }

    for (k=0; (k<fak[candidates[0][0]]) && !(*abbruch); k++)
    { 
    if (k)
    { zw1=candidates[0][vert[k][0]];
    zw2=candidates[0][vert[k][1]];
    puffer=bild[zw1];
    bild[zw1] = bild[zw2];
    bild[zw2]=puffer;
    if (suche==zw1) suche=zw2; else { if (suche==zw2) suche=zw1; } 
                 /* hat sich der Wert des naechsten Knotens veraendert ? */
    }
    normalweiter(graph,adj,bild,suche,naechstenummer,abbruch); 
    } /* Ende for [0][0]*/
    } /* Ende for [1][0]*/
    } /* Ende for [2][0]*/

   } /* Ende if vertex<.....*/

    for (j=1; j<=candidates[2][0]; j++) bild[candidates[2][j]]=unbelegt;
    for (j=1; j<=candidates[1][0]; j++) bild[candidates[1][j]]=unbelegt;
    for (j=1; j<=candidates[0][0]; j++) bild[candidates[0][j]]=unbelegt;
    for (j=1; j<=einserspeicher[0]; j++) bild[einserspeicher[j]]=unbelegt;

   }/* Ende es gibt ein Loch */
   else /* d.h. kein loch */
     { 
     co=compare(bild,graph[wo],graph[bild[wo]]);
     if (co<0) *abbruch=1;
     else
     if ((vertex<graph[0][0])&&(co==0))
       {
        for (suche=1; bild[suche]!=vertex+1; suche++); 
                               /*sucht den Knoten mit bild vertex+1 */
        normalweiter(graph,adj,bild,suche,naechstenummer,abbruch); }
     }
 } /* ende not drin[vertex] */
}/* ende adj[wo]>0 denn wenn die adjazenz 0 ist, muss nichts gemacht werden */

} /* ende normalweiter */


/******************************NAECHSTERSUCHE***********************/

int naechstersuche(GRAPH graph, ADJAZENZ adj, unsigned char wo,
                   BILD bild, unsigned char *naechster)

//GRAPH graph;
//ADJAZENZ adj;
//unsigned char wo; /* der wirkliche name des zu behandelnden knotens */
//BILD bild;
//unsigned char *naechster; /* die adresse fuer den naechsten knoten, an dem
//                             wahlmoeglichkeit besteht */
/* sucht den naechsten zu bearbeitenden Knoten, startend bei wo. Der naechste
   Knoten ist dadurch gekennzeichnet, dass er oder sein "Original" noch nicht
   fertig sind, oder dass in seinem Bild noch unbelegte Felder sind. Der 
   Returnwert ist der Vergleich bis zu diesem Knoten -- aber:
   Ist der letzte knoten an dem noch wahlmoeglichkeit besteht noch nicht
   voll saturiert, so wird nur ein negatives testergebnis mit einbezogen
   in die auswertung von test -- ein positives kann sich noch aendern    */
{
int test,puf,echt,ende,j;

ende=test=0;
echt = bild[wo];
for (j=0;(j<adj[wo]) && !ende;j++) 
                          if (bild[graph[wo][j]]==unbelegt) ende=1;

for ( ; fertig[wo] && !ende && fertig[echt] && !test && (echt<graph[0][0]); )
           { 
             test = compare(bild,graph[wo],graph[echt]);
             echt++;
             for (wo=1; bild[wo]!=echt; wo++);
             for (j=0;(j<adj[wo]) && !ende;j++) 
                          if (bild[graph[wo][j]]==unbelegt) ende=1;
           }
/* Die wo-echt Kombination hat jetzt den ersten wert wo etwas schief geht
   und test den Vergleichswert davor */

if (test) return(test); /* dieses Testergebnis gilt immer */
                        /* *naechster ist aber nicht sinnvoll belegt */
    else
       if (!fertig[wo] || ende || !fertig[echt]) 
                         { *naechster=wo;
                           puf= compare(bild,graph[wo],graph[echt]);
                           if (puf<0) return(puf); else return(0);
                         }
       else /* d.h. echt == graph[0][0] */
	   { *naechster=wo; return(0); }

}


/**************************TOTSETZEN**********************************/

void totsetzen(BRANCH *zweig, int woher, unsigned char tiefe)

//BRANCH *zweig;  /* welcher zweig soll bearbeitet werden */
//int woher;  /* woher kommt man ? */
//unsigned char tiefe;  /* das zeichen, ab wann ein zweig tot ist */

/* diese funktion laeuft den baum rueckwaerts und setzt den totzeiger in
   die richtung aus der man kommt auf tiefe. Wenn der andere zeiger auch
   blockiert ist, geht es noch einen zweig zurueck, ...                */

{

(zweig->tot)[woher]=tiefe;
(zweig->anzahltot)++;
if (zweig->anzahltot == zweig->anzahlnachfolger)
  { if (zweig->vorgaenger != nil) 
             totsetzen(zweig->vorgaenger,zweig->zweignummer,tiefe);
    else (roots[zweig->rootnumber]).wurzeltot=tiefe; }
}


/*************************ZWEIGINITIALISIEREN***************************/

void zweiginitialisieren(BRANCH *zweig, BRANCH *sollvorgaenger, BILD sollbild,
                         unsigned char sollnextbild, unsigned char solltiefe,
                         unsigned char sollwo, int sollzweignummer,
                         unsigned char rootnb)

/* Beschreibt ein Zweigelement als Blatt. Nachfolger werden nicht
   eingefuegt ( das passiert in Baumbau ) */

//BRANCH *zweig; /* der zu initialisierende zweig */
//BRANCH *sollvorgaenger;
//BILD sollbild;
//unsigned char sollnextbild;
//unsigned char solltiefe;
//unsigned char sollwo; 
//int sollzweignummer;
//unsigned char rootnb;
{

zweig->blatt = 1;
memcpy(zweig->bild,sollbild,sizeof(BILD));
zweig->vorgaenger=sollvorgaenger;
zweig->nextbild=sollnextbild;
zweig->tiefe = solltiefe;
zweig->zweignummer=sollzweignummer;
zweig->rootnumber=rootnb;
zweig->anzahltot=zweig->anzahlnachfolger=0;
zweig->tot=nil;
zweig->successors=nil;
zweig->wo = sollwo;


}/* ende zweiginitialisieren */

/***********************************BAUMBAU*****************************/


void baumbau(GRAPH graph, ADJAZENZ adj, BRANCH *zweig, int *abbruch)

//GRAPH graph;
//ADJAZENZ adj;
//BRANCH *zweig;
//int *abbruch;

/* baumbau versucht, startend bei zweig den baum auszubauen -- 
   zweig selbst bleibt dabei unveraendert bis auf die Zeigerwerte */

{
int i,j,k,suche,co,zaehler;
unsigned char merke, naechstenummer, wo, sollwo;
int zw1, zw2, puffer, puf;
unsigned char *(candidates[3]);
unsigned char kandidatenspeicher[9][9];
unsigned char einserspeicher[9];
int maxmultiplicity, numberofcand;
BILDELEMENT *bild;

wo=zweig->wo;
bild= zweig->bild;
if (fertig[wo] && fertig[bild[wo]])
{
(blaetterzahl[zweig->rootnumber])--; 
zweig->blatt=0;  /* auf jeden Fall kein Blatt mehr */

einserspeicher[0]=0;
naechstenummer=zweig->nextbild;

   maxmultiplicity=0;
   for (i=0; i<9; i++) kandidatenspeicher[i][0]=0;
   suche= -1;
   for (i=0; i<adj[wo]; )
      { if (bild[graph[wo][i]]==unbelegt)
        { if ((adj[graph[wo][i]]==1) && (fertig[graph[wo][i]]))
           {
            (einserspeicher[0])++;
            einserspeicher[einserspeicher[0]]=graph[wo][i];
            i++;
            suche=1;
           }
         else
          { suche=i; merke=graph[wo][i]; i++;
            for( ;(graph[wo][i]==merke) && (i<adj[wo]); i++); 
            suche=i-suche; kandidatenspeicher[suche][0]++; 
            kandidatenspeicher[suche][kandidatenspeicher[suche][0]]=merke; 
            if (suche>maxmultiplicity) maxmultiplicity=suche; }
	}
         else i++; }
     
if (suche>=0) /* d.h. es gibt ein Loch */
  {
     numberofcand=0;
      candidates[0]=candidates[1]=candidates[2]=kandidatenspeicher[0]; 
         /* Diese werte werden nicht notwendig belegt, muessen aber fuer die 
            Abfrage in der for-schleife den Defaultwert 0 beim Element c[i][0]
            haben */
    for (i=maxmultiplicity; i>0; i--)
      { if (kandidatenspeicher[i][0]) 
                            { candidates[numberofcand]=kandidatenspeicher[i];
                              numberofcand++;
                              for (j=1; j<=kandidatenspeicher[i][0]; j++) 
                                 {
                                 bild[kandidatenspeicher[i][j]]=naechstenummer;
                                 naechstenummer++; }
                                      }
      }

for (j=1; j<=einserspeicher[0]; j++)
    { bild[einserspeicher[j]]=naechstenummer;
      naechstenummer++;
    }


    co=compare(bild,graph[wo],graph[bild[wo]]);
    if (co<0) *abbruch=1;
    else 
    if (co>0) { totsetzen(zweig->vorgaenger,zweig->zweignummer,graph[0][1]);}

    else /* d.h. co==0 */
    {
    zaehler=0;
    zweig->anzahlnachfolger = fak[candidates[0][0]]*fak[candidates[1][0]]
                              *fak[candidates[2][0]];
    zweig->successors = (struct zweig **)calloc(zweig->anzahlnachfolger,
                                         sizeof(struct zweig *));
    zweig->tot = (unsigned char *)calloc(zweig->anzahlnachfolger,
                                         sizeof(unsigned char));
    puf=bild[wo]+1;
    for (suche=1; bild[suche]!=puf; suche++); /*das urbild von bild[wo]+1*/

    for (i=0; (i<fak[candidates[2][0]]) && !(*abbruch); i++)
    { 
    if (i)
    { zw1=candidates[2][vert[i][0]];
    zw2=candidates[2][vert[i][1]];
    puffer=bild[zw1];
    bild[zw1] = bild[zw2];
    bild[zw2]=puffer;
    if (suche==zw1) suche=zw2; else { if (suche==zw2) suche=zw1; } 
             /* hat sich der Wert des naechsten Knotens veraendert ? */
    }

    for (j=0; (j<fak[candidates[1][0]]) && !(*abbruch); j++)
    { 
    if (j)
    { zw1=candidates[1][vert[j][0]];
    zw2=candidates[1][vert[j][1]];
    puffer=bild[zw1];
    bild[zw1] = bild[zw2];
    bild[zw2]=puffer;
    if (suche==zw1) suche=zw2; else { if (suche==zw2) suche=zw1; } 
             /* hat sich der Wert des naechsten Knotens veraendert ? */
    }

    for (k=0; (k<fak[candidates[0][0]]) && !(*abbruch); k++)
    { 
    if (k)
    { zw1=candidates[0][vert[k][0]];
    zw2=candidates[0][vert[k][1]];
    puffer=bild[zw1];
    bild[zw1] = bild[zw2];
    bild[zw2]=puffer;
    if (suche==zw1) suche=zw2; else { if (suche==zw2) suche=zw1; } 
             /* hat sich der Wert des naechsten Knotens veraendert ? */
    }
    puf = naechstersuche(graph,adj,suche,bild,&sollwo);

    if (puf<0) *abbruch=1;
    else
    if (puf>0) {
                 totsetzen(zweig,zaehler,graph[0][1]);
                 zaehler++;
               }

    else /* d.h. puf==0 */
    {
    (blaetterzahl[zweig->rootnumber])++;
    zweig->successors[zaehler]=
                (struct zweig *)malloc(sizeof(struct zweig));
    zweiginitialisieren(zweig->successors[zaehler],zweig,bild,
                  naechstenummer,graph[0][1],sollwo,zaehler,zweig->rootnumber);
    baumbau(graph,adj,zweig->successors[zaehler],abbruch);
    zaehler++;
    }
    } /* Ende for [0][0]*/
    } /* Ende for [1][0]*/
    } /* Ende for [2][0]*/

    } /* Ende if vertex<.....*/
    for (j=1; j<=candidates[2][0]; j++) bild[candidates[2][j]]=unbelegt;
    for (j=1; j<=candidates[1][0]; j++) bild[candidates[1][j]]=unbelegt;
    for (j=1; j<=candidates[0][0]; j++) bild[candidates[0][j]]=unbelegt;
    for (j=1; j<=einserspeicher[0]; j++) bild[einserspeicher[j]]=unbelegt;

  } /* ende if es gibt loecher */
  else /* d.h. keine loecher */
    {
    puf = naechstersuche(graph,adj,wo,bild,&sollwo);
    if (puf<0) *abbruch=1;
    else
    if (puf>0) totsetzen(zweig->vorgaenger,zweig->zweignummer,graph[0][1]);
    else /* d.h. puf==0 */
    if (bild[wo]<graph[0][0]) /* ansonsten ist man fertig */
    {
    (blaetterzahl[zweig->rootnumber])++;

    zweig->anzahlnachfolger=1;
    zweig->successors= (struct zweig **)malloc(sizeof(struct zweig*));
    zweig->tot=(unsigned char *)calloc(1,sizeof(unsigned char));
    zweig->successors[0]=
                (struct zweig *)malloc(sizeof(struct zweig));
    zweiginitialisieren(zweig->successors[0],zweig,bild,
                        naechstenummer,graph[0][1],sollwo,0,zweig->rootnumber);
    baumbau(graph,adj,zweig->successors[0],abbruch);
    }
    else { zweig->blatt=1; (blaetterzahl[zweig->rootnumber])++; }
    } /*ende else -- d.h. keine loecher */
  } /* ende fertig[wo] && fertig[bild[wo]] */
  else
   normalweiter(graph,adj,zweig->bild,zweig->wo,zweig->nextbild,abbruch);

}/* ende baumbau */




/***************************ENDMINI*****************************************/


int endmini (GRAPH graph, ADJAZENZ adj)

//GRAPH graph;
/* graph[0][0] ist die anzahl der knoten, graph[0][1] die anzahl der kanten
   ( also tiefe )  */

//ADJAZENZ adj; /* adj[i] ist der grad von knoten i */


{
int i,j,abbruch;
BILD image;
BRANCH *zweig;

if (orig_non[1]>=7) preprocessing(graph,adj); /* nur wenn es sich lohnt */


for (i=1; i<=graph[0][0]; i++) image[i]=unbelegt;

abbruch = 0;

for (i=1; (i<=graph[0][0]) && (!abbruch); i++)
             { if (roots[i].tiefe) /*d.h.: schon initialisiert*/
                  {
                     for (j=0; (j<blaetterzahl[i]) && (!abbruch); j++)
		       { zweig=blaetterkette[i][j]; 
                         normalweiter(graph,adj,zweig->bild,zweig->wo,
                                              zweig->nextbild,&abbruch); 
		       }
		  }
                        else if (adj[i]>1)
                         {
                  image[i]=1;
                  normalweiter(graph,adj,image,i,2,&abbruch);
                  image[i]=unbelegt;
		         } /* Ende else */


} /* Ende grosse for-Schleife */

/*fprintf(stderr,"Ergebnis in endmini: %d\n",!abbruch);*/

return (!abbruch);
}
    

/**************************EINFUGEN************************************/

void einfugen (GRAPH graph, ADJAZENZ adj, unsigned char v, unsigned char w)
//GRAPH graph;
//ADJAZENZ adj;
//unsigned char v, w;
/* Fuegt die Kante (v,w) in den Graphen graph ein. Dabei wird aber davon */
/* ausgegangen, dass in adj die wirklich aktuellen werte fuer die */
/* Adjazenzen stehen. Die adjazenzen werden dann aktualisiert. */

{ 
graph[v][adj[v]]=w;
graph[w][adj[w]]=v;
adj[v]++;
adj[w]++;
stillness[adj[w]]--;
stillness[adj[v]]--;
graph[0][1]++;
}

/*************************ENTFERNEN*************************************/

void entfernen (GRAPH graph, ADJAZENZ adj, unsigned char v, unsigned char w)
//GRAPH graph;
//ADJAZENZ adj;
//unsigned char v, w;
/* Entfernt die Kante (v,w) aus dem Graphen graph. Dabei wird aber davon */
/* ausgegangen, dass in adj die wirklich aktuellen werte fuer die */
/* Adjazenzen stehen. Die adjazenzen werden dann aktualisiert. */

{ 
graph[v][adj[v]-1]=leer;
graph[w][adj[w]-1]=leer;
stillness[adj[w]]++;
stillness[adj[v]]++;
adj[v]--;
adj[w]--;
graph[0][1]--;
}




/******************************KETTENENTFERNEN***********************/


void kettenentfernen(unsigned char bis_wo )

//unsigned char bis_wo;

{
int i;

for (i=1; i<= bis_wo; i++)
      if (blaetterkette[i] != nil) { free(blaetterkette[i]);
                                     blaetterkette[i]=nil;
				   }
}

/**********************************SUCHE*****************************/

void suche(GRAPH graph, ADJAZENZ adj, BRANCH *zweig, int *abbruch)

//GRAPH graph;
//ADJAZENZ adj;
//BRANCH *zweig;
//int *abbruch;

/* suche durchlaeuft einfach den Baum auf der suche nach Blaettern */


{ int i;

if (zweig->blatt) baumbau(graph,adj,zweig,abbruch);

else
     for (i=0; i<zweig->anzahlnachfolger; i++)
     if ((!zweig->tot[i]) && !(*abbruch))
               suche(graph,adj,zweig->successors[i],abbruch); 
}

/*****************************ROOTINITIALISIEREN************************/


void rootinitialisieren(GRAPH graph, ADJAZENZ adj, ROOT *which,
                        int number, int *abbruch)

//GRAPH graph;
//ADJAZENZ adj;
//ROOT *which;  /* zeiger auf die zu bearbeitende wurzel */
//int number; /* welcher knoten soll die 1 als bild bekommen */
//int * abbruch;

/* Die Struktur dieser Funktion ist sehr aehnlich wie die von
   normalweiter und baumbau. Die jeweiligen Unterschiede sind
   ( besonders zwischen rootinitialisieren und baumbau ) gering */

{
int i,j,k,suche,co,zaehler;
unsigned char merke, naechstenummer, sollwo;
int zw1, zw2, puffer;
unsigned char *(candidates[3]); /* Die Loecherkandidaten, fuer die kleinste 
        nummer werden in dem vector, auf den candidates[0],... zeigt, 
        gespeichert und entsprechend weiter. Bei Auslegung auf Maximalvalenz 
        8 reichen 3 Zeiger -- ueberlegen*/
unsigned char kandidatenspeicher[9][9]; /* In Kandidatenspeicher[i][0] wird 
        die Anzahl der Kandidaten mit multiplizitaet[i] gespeichert und in 
        den folgenden elementen die kandidaten selbst */
unsigned char einserspeicher[9];
int maxmultiplicity, numberofcand;
BILDELEMENT *bild;
struct zweig * zweige;

zweige= &(which->zweige);
bild= zweige->bild;
for (i=1; i<=knoten; i++) bild[i]=unbelegt;
bild[number]=1;
zweige->wo=number;
zweige->tiefe=graph[0][1]; /* Die beiden vorigen Zuweisungen sind an sich
                              ueberfluessig und dienen nur der Uebersicht-
                              lichkeit des Baumes */
naechstenummer=2;
blaetterzahl[number]=0;
which->tiefe=graph[0][1];
zweige->vorgaenger=nil;
zweige->successors=nil;
zweige->tot=nil;
zweige->rootnumber=number;
einserspeicher[0]=0;

   maxmultiplicity=0;
   for (i=0; i<9; i++) kandidatenspeicher[i][0]=0;
   suche= -1;
   for (i=0; i<adj[number]; )
      { if (bild[graph[number][i]]==unbelegt)
        { if ((adj[graph[number][i]]==1) && (fertig[graph[number][i]]))
           {
            (einserspeicher[0])++;
            einserspeicher[einserspeicher[0]]=graph[number][i];
            i++;
            suche=1;
           }
         else
          { suche=i; merke=graph[number][i]; i++;
            for( ;(graph[number][i]==merke) && (i<adj[number]); i++); 
            suche=i-suche; kandidatenspeicher[suche][0]++; 
            kandidatenspeicher[suche][kandidatenspeicher[suche][0]]=merke; 
            if (suche>maxmultiplicity) maxmultiplicity=suche; }
	  }
         else i++; }

/* Hier gibt es immer mindestens ein loch */
     numberofcand=0;
      candidates[0]=candidates[1]=candidates[2]=kandidatenspeicher[0]; 
            /* Diese werte werden nicht notwendig belegt, muessen aber fuer 
            die Abfrage in der for-schleife den Defaultwert 0 beim Element
            c[i][0] haben */
    for (i=maxmultiplicity; i>0; i--)
      { if (kandidatenspeicher[i][0]) 
                            { candidates[numberofcand]=kandidatenspeicher[i];
                              numberofcand++;
                              for (j=1; j<=kandidatenspeicher[i][0]; j++) 
                                 {
                                 bild[kandidatenspeicher[i][j]]=naechstenummer;
                                 naechstenummer++; }
                                      }
      }

for (j=1; j<=einserspeicher[0]; j++)
    { bild[einserspeicher[j]]=naechstenummer;
      naechstenummer++;
    }


    co=compare(bild,graph[number],graph[1]);

    if (co<0) *abbruch=1;
    else 
    if (co>0) which->wurzeltot = graph[0][1];
         /* which muss ja ein fertiger Knoten sein und von daher aendert
            sich dies Verhaeltnis nicht mehr */
    else
    {
    which->wurzeltot=zweige->anzahltot=0;
    zaehler=0;

    zweige->anzahlnachfolger = fak[candidates[0][0]]*fak[candidates[1][0]]
                              *fak[candidates[2][0]];

    zweige->successors = (struct zweig **)calloc(zweige->anzahlnachfolger,
                                         sizeof(struct zweig *));
    zweige->tot = (unsigned char *)calloc(zweige->anzahlnachfolger,
                                         sizeof(unsigned char));
    suche=candidates[0][1]; /* das urbild von 2 */

    for (i=0; (i<fak[candidates[2][0]]) && !(*abbruch); i++)
    { 
    if (i)
    { zw1=candidates[2][vert[i][0]];
    zw2=candidates[2][vert[i][1]];
    puffer=bild[zw1];
    bild[zw1] = bild[zw2];
    bild[zw2]=puffer;
    if (suche==zw1) suche=zw2; else { if (suche==zw2) suche=zw1; } 
             /* hat sich der Wert des naechsten Knotens veraendert ? */
    }

    for (j=0; (j<fak[candidates[1][0]]) && !(*abbruch); j++)
    { 
    if (j)
    { zw1=candidates[1][vert[j][0]];
    zw2=candidates[1][vert[j][1]];
    puffer=bild[zw1];
    bild[zw1] = bild[zw2];
    bild[zw2]=puffer;
    if (suche==zw1) suche=zw2; else { if (suche==zw2) suche=zw1; } 
             /* hat sich der Wert des naechsten Knotens veraendert ? */
    }

    for (k=0; (k<fak[candidates[0][0]]) && !(*abbruch); k++)
    { 
    if (k)
    { zw1=candidates[0][vert[k][0]];
    zw2=candidates[0][vert[k][1]];
    puffer=bild[zw1];
    bild[zw1] = bild[zw2];
    bild[zw2]=puffer;
    if (suche==zw1) suche=zw2; else { if (suche==zw2) suche=zw1; } 
             /* hat sich der Wert des naechsten Knotens veraendert ? */
    }
    puffer = naechstersuche(graph,adj,suche,bild,&sollwo);
    if (puffer<0) *abbruch=1;
    else
    if (puffer>0) { zweige->tot[zaehler]=graph[0][1];
                    zweige->anzahltot++;
                    zaehler++;
                  }
    else /* d.h. puffer==0 */
    {
    blaetterzahl[number]++;
    zweige->successors[zaehler]=
                (struct zweig *)malloc(sizeof(struct zweig));
    zweiginitialisieren(zweige->successors[zaehler],zweige,bild,
                            naechstenummer,graph[0][1],sollwo,zaehler,number);
    baumbau(graph,adj,zweige->successors[zaehler],abbruch);
    zaehler++;
    }
    } /* Ende for [0][0]*/
    } /* Ende for [1][0]*/
    } /* Ende for [2][0]*/

    if (zweige->anzahltot==zweige->anzahlnachfolger) 
                                    which->wurzeltot=graph[0][1];

   } /* Ende if vertex<.....*/
    for (j=1; j<=candidates[2][0]; j++) bild[candidates[2][j]]=unbelegt;
    for (j=1; j<=candidates[1][0]; j++) bild[candidates[1][j]]=unbelegt;
    for (j=1; j<=candidates[0][0]; j++) bild[candidates[0][j]]=unbelegt;
    for (j=1; j<=einserspeicher[0]; j++) bild[einserspeicher[j]]=unbelegt;
}


/*************************MINIDAR************************************/

int minidar (GRAPH graph, ADJAZENZ adj, unsigned char vertex )

//GRAPH graph;
/* graph[0][0] ist die anzahl der knoten, graph[0][1] die anzahl der kanten
   ( also tiefe ) */

//ADJAZENZ adj; /* adj[i] ist der grad von knoten i */
//unsigned char vertex;

{
int i,abbruch;
BILD image;

if (orig_non[1]>=7) preprocessing(graph,adj); /* nur wenn es sich lohnt */
/* Achtung ! Nach einmaligem Aufruf von preprocessing sind die drin[]-Werte
   korrupt und werden erst beim naechsten Aufruf wieder in Ordnung gebracht */

for (i=1; i<=graph[0][0]; i++) 
      { fertig[i]= ( !((drin[i][0]) | (drin[i][1])) && 
                    ((i<vertex) || (stillness[adj[i]+1]==0)));
                                 image[i]=unbelegt; }

abbruch = 0;

for (i=1; (i<=graph[0][0]) && (!abbruch); i++)
             if (fertig[i])  /* nur fertige Knoten koennen als Wurzel dienen */
             { if (roots[i].tiefe) /*d.h.: schon initialisiert*/
	       { if (!(roots[i].wurzeltot))
                    suche(graph,adj,&(roots[i].zweige),&abbruch); }
                        else /* d.h. root noch nicht initialisiert */
                          rootinitialisieren(graph,adj,roots+i,i,&abbruch);
                          /* rootinitialisieren ruft auch baumbau auf */
	     }
             else /* d.h. der Knoten ist noch nicht fertig */
                         {
                  image[i]=1;
                  normalweiter(graph,adj,image,i,2,&abbruch);
                  image[i]=unbelegt;
		         } /* Ende else */

/*fprintf(stderr,"Ergebnis in minidar: %d\n",!abbruch);*/

return (!abbruch);
}
    

/**************************MEMOCOMP*************************************/


int memocomp(char *adr1, char *adr2, int anzahl)
/* Der Ersatz fuer das auf einigen dec-rechnern nicht funktionierende memcmp */

//unsigned char *adr1, *adr2;
//int anzahl;

{
for ( ;(*adr1 == *adr2) && (anzahl>1); adr1++, adr2++, anzahl--);
return ((int)*adr1 - (int)*adr2);
/* geaendert 11.2.2004 */
}


/*****************************KONSTRUKT**********************************/



void konstrukt(GRAPH graph, ADJAZENZ adj, unsigned char vertex)
/* konstruiert weiter am graphen und zwar am Knoten vertex.*/

//GRAPH graph;
//ADJAZENZ adj;
//unsigned char vertex;

{
long int anfang, i, test,baumgebaut, kantenvorher;
unsigned char analysed; /* wurde eine baumanalyse gemacht ? */

/*dbz++; fprintf(stderr," dbz: %d   vertex: %d\n",dbz,vertex);
fprintf(stderr,"Anfang konstrukt\n");
schreibegraph(graph);*/


if (mod && (graph[0][1]==splitlevel))
    { splitzaehler++;
      if (splitzaehler==mod) splitzaehler=0;
      if (splitzaehler != rest) return; }


if ((recover) && (recover==vertex) &&
        (!memocomp((unsigned char *)graph,(unsigned char *)lastgraph, sizeof(GRAPH))) )
       /* d.h. der lastgraph wurde rekonstruiert */
                { recover=0; }
else

{

analysed=0;

kantenvorher=graph[0][1];
if ((kantenvorher==treeborder+1) && (blaetterzahl[0]==0))
    { blaetterzahl[0]=kantenvorher; 
      analysed=1;
      baumanalyse(graph[0][0]); }

if (vertex>knotenzahl) /* d.h.: Der graph ist fertig */
                  { 
		    if (endmini(graph,adj)) aufschreiben(graph,adj);
		  } 
else
   /* (adj[vertex]==0) wuerde bedeuten, dass der Rest nicht mit dem bereits
      konstruierten zusammenhaengen kann */
    
   if (adj[vertex]>0)
     {  
      /* Zunaechst wird versucht, einfach weiter zu machen: */
      if (nr_of_nodes[adj[vertex]]) 
                /* d.h.: Dieser Knoten hat eine "nuetzliche" Valenz */
            if ( (!recover) || (lastgraph[vertex][adj[vertex]]==leer) )
                { if (mod && (graph[0][1]==splitlevel)) splitzaehler--;
		  /* ansonsten wuerde splitzaehler hochgezaehlt, ohne eine Kante eingefuegt zu
		     haben, was dazu fuehren wuerde, dass (mindestens) zwei verschiedene
		     Werte fuer splitzaehler vorhanden waeren -- und bei einem wuerde er mit
		     Sicherheit IN JEDEM FALL verworfen -- also Fehler. */
                  nr_of_nodes[adj[vertex]]--;
                  konstrukt(graph,adj,vertex+1);
                  nr_of_nodes[adj[vertex]]++;
                }

      if (stillness[adj[vertex]+1]) /* es werden noch Knoten mit groesserer */
                                    /* Valenz gebraucht                     */
            {
            anfang=graph[vertex][adj[vertex]-1]; 
	                                 /*geaendert fuer Multigraphen*/ 
            if (anfang<vertex) anfang=vertex+1;
                   if ((adj[vertex]>=multiplicity) && 
                            (graph[vertex][adj[vertex]-multiplicity]==anfang))
                                 anfang++;
          /* sonst haette dieser Knoten eine groessere multiplizitaet als 1 */
                   if (zykelkanten<maxzykelkanten)
                       for (i=anfang; (i<=graph[0][0]); i++)
                           /* Hier wird jedesmal ein Zykel geschlossen ! */
		      { if (adj[vertex]==adj[i]) test=(stillness[adj[i]+1]>=2);
                                   else test=(stillness[adj[i]+1]>=1);
                           /* stillness[adj[vertex]+1] ist hier ja eh true */
			if (bipartite) (test=test && (colour[vertex] != colour[i]));
                  test=test && abstandok(graph,vertex,(unsigned char)i,taille);
                        if ( test )
                               {
                                 einfugen(graph,adj,vertex,i);
                                 zykelkanten++;
                                 baumgebaut=0;
                                 if (recover) test=
                                         (memocomp(graph[1],lastgraph[1],
                                        (knoten+1)*(vertex-1)+adj[vertex])==0);
                                                 else test=1;
/* Achtung ! Die vorige Zeile benutzt explizit eine vorgegebene Form von
   GRAPH -- bei einer anderen Definition kommt es zu Fehlern ! */

                                   if (test)
				       { if (kantenvorher<=treeborder)
					     { baumgebaut = 1;
                                              test = minidar(graph,adj,vertex);
					     }
                                         else 
                                     if (kantenvorher<=border)
                                                  test=endmini(graph,adj);
				        }

                                 if (test) konstrukt(graph,adj,vertex);
                                 if (baumgebaut)
                                 baumaufraeumen(graph[0][1],graph[0][0]);
                                 entfernen(graph,adj,vertex,i);
                                 zykelkanten--;
                               }
      	       } /* ende for "alle bereits enthaltenen groesseren Knoten */
       /* Jetzt noch neue Knoten hinzufuegen */
                   if (graph[0][0]<knotenzahl)
                       {  graph[0][0]++; 
			  if (bipartite) { if (colour[vertex]==1) 
					     { colour[graph[0][0]]=0;
					       (anzcolours[0])++; }
					   else 
					     { colour[graph[0][0]]=1;
					       (anzcolours[1])++; }
					 }
                          einfugen(graph,adj,vertex,graph[0][0]);
                          baumgebaut=0;
                          if (recover) test=(memocomp(graph[1],lastgraph[1],
                                        (knoten+1)*(vertex-1)+adj[vertex])==0);
                                                 else test=1;
                                             /* siehe oben */
			  if (bipartite && regular)
			    { if ((anzcolours[0]>knotenzahl/2) || (anzcolours[1]>knotenzahl/2)) test=0; }
			  
			  if (test)
			    { if (kantenvorher<=treeborder)
				{ baumgebaut = 1;
				  test = minidar(graph,adj,vertex);
				}
			    else 
			      if (kantenvorher<=border)
				test=endmini(graph,adj);
			    }
			  
                          if (test) konstrukt(graph,adj,vertex);
                          if (baumgebaut)
                          baumaufraeumen(graph[0][1],graph[0][0]);
                          entfernen(graph,adj,vertex,graph[0][0]);
			  (anzcolours[colour[graph[0][0]]])--;
                          graph[0][0]--;
                        }
               } /* ende if stillness */

if (analysed) { kettenentfernen(graph[0][0]); blaetterzahl[0]=0; }

} /* ende if adjazenz[vertex]>0 */

/*if (!recover)
  {  
     if ( time(0) >= Sicherungszeit )
       { wegspeichern(globalliste,graphenzahl,graph,vertex); graphenzahl=0;
         Sicherungszeit=time(0)+Sicherungsintervall; }
  }*/
/* der time(0) aufruf braucht immens viel Systemzeit wenn er, wie hier,
   oft gemacht wird -- auf signale wie in minibaum umstellen */


} /* ende else nach lastgraph wiedergefunden */


} /* ende konstruct */





/*****************************V1_KONSTRUKT**********************************/



void v1_konstrukt(GRAPH graph, ADJAZENZ adj, unsigned char vertex)
/* Beginnt die Konstruktion bei vertex 1. Im wesentlichen identisch mit
   konstrukt -- hier wird lediglich z.B. maxmultiplicity festgelegt    */

//GRAPH graph;
//ADJAZENZ adj;
//unsigned char vertex;

{
long int anfang, i, test,baumgebaut, kantenvorher, altmult;
unsigned char analysed;

/* Hier wird recover nie auf 0 gesetzt, da hier auch nie Sicherungen
   gemacht werden */

analysed=0;

kantenvorher=graph[0][1];
if ((kantenvorher==treeborder+1) && (blaetterzahl[0]==0))
    { blaetterzahl[0]=kantenvorher;
      analysed=1;
      baumanalyse(graph[0][0]); }

      /* Zunaechst wird versucht, einfach weiter zu machen: */
      if (nr_of_nodes[adj[vertex]]) 
                /* d.h.: Dieser Knoten hat eine "nuetzliche" Valenz */
            if ( (!recover) || (lastgraph[vertex][adj[vertex]]==leer) )
                { 
                  nr_of_nodes[adj[vertex]]--;
                  konstrukt(graph,adj,vertex+1);
                  nr_of_nodes[adj[vertex]]++;
                }

      if (stillness[adj[vertex]+1]) /* es werden noch Knoten mit groesserer */
                                    /* Valenz gebraucht                     */
            {
            anfang=graph[vertex][adj[vertex]-1]; 
	                                 /*geaendert fuer Multigraphen*/ 
            if (taille>=3) anfang++;
            if (anfang<=vertex) anfang=vertex+1;
                   if (zykelkanten<maxzykelkanten)
                       for (i=anfang; (i<=graph[0][0]); i++)
                           /* Hier wird jedesmal ein Zykel geschlossen !
                              In diesem Fall eine Doppelkante */
		      { if (adj[vertex]==adj[i]) test=(stillness[adj[i]+1]>=2);
                                   else test=(stillness[adj[i]+1]>=1);
                           /* stillness[adj[vertex+1]] ist hier ja eh true */
                               if (i>2) 
                           test = (test && (graph[1][adj[1]-multiplicity]!=i));
			   if (bipartite) (test=test && (colour[vertex] != colour[i]));
                           if (i==2) test = test && ( (adj[1] < maxvalence-1) 
                                                      || (knotenzahl==2) );
                  /* keine multiplizitaet darf hoeher sein, als die von 1->2 */
                        if ( test )
                               { 
                                 einfugen(graph,adj,vertex,i);
                                 altmult=multiplicity;
                                 if (i==2) multiplicity=adj[1];
                                 zykelkanten++;
                                 baumgebaut=0;
                                 if (recover) test=
                                         (memocomp(graph[1],lastgraph[1],
                                        (knoten+1)*(vertex-1)+adj[vertex])==0);
                                                 else test=1;
/* Achtung ! Die vorige Zeile benutzt explizit eine vorgegebene Form von
   GRAPH -- bei einer anderen Definition kommt es zu Fehlern ! */

                                   if (test)
				       { if (kantenvorher<=treeborder)
					     { baumgebaut = 1;
                                              test = minidar(graph,adj,vertex);
					     }
                                         else 
                                      if (kantenvorher<=border)
                                                  test=endmini(graph,adj);
				        }

                                 if (test) v1_konstrukt(graph,adj,vertex);
			                   
                                 if (baumgebaut)
                                 baumaufraeumen(graph[0][1],graph[0][0]);
                                 entfernen(graph,adj,vertex,i);
                                 multiplicity=altmult;
                                 zykelkanten--;
                               }
		       } 
                   /* ende for "alle bereits enthaltenen groesseren Knoten */
       /* Jetzt noch neue Knoten hinzufuegen */ 
                   if (graph[0][0]<knotenzahl)
                       {  graph[0][0]++; 
			  if (bipartite) { if (colour[vertex]==1) 
					     { colour[graph[0][0]]=0;
					       (anzcolours[0])++; }
					   else 
					     { colour[graph[0][0]]=1;
					       (anzcolours[1])++; }
					 }
                          einfugen(graph,adj,vertex,graph[0][0]);
                          baumgebaut=0;
                          if (recover) test=(memocomp(graph[1],lastgraph[1],
                                        (knoten+1)*(vertex-1)+adj[vertex])==0);
                                                 else test=1;
                                             /* siehe oben */
			  if (bipartite && regular)
			    { if ((anzcolours[0]>knotenzahl/2) || (anzcolours[1]>knotenzahl/2)) test=0; }
			  if (test)
			    { if (kantenvorher<=treeborder)
				{ baumgebaut = 1;
				  test = minidar(graph,adj,vertex);
				}
			    else 
                                      if (kantenvorher<=border)
					test=endmini(graph,adj);
			    }
			  
                          if (test) v1_konstrukt(graph,adj,vertex);
                          if (baumgebaut)
			    baumaufraeumen(graph[0][1],graph[0][0]);
                          entfernen(graph,adj,vertex,graph[0][0]);
			  (anzcolours[colour[graph[0][0]]])--;
                          graph[0][0]--;
                        }

               } /* ende if stillness */

if (analysed) {  kettenentfernen(graph[0][0]); blaetterzahl[0]=0; }

} /* ende v1_konstrukt */



/**************************INITIALIZE*********************************/

void initialize(GRAPH graph, ADJAZENZ adj)
//GRAPH graph;
//ADJAZENZ adj;

/* initialisiert graph, adj und roots*/
/* Achtung ! hier wird die genaue definition von GRAPH benutzt !!!!!!!*/

{
int i,j;

blaetterzahl[0]=adj[0]=0;

for (i=1; i<=knoten; i++) { for (j=0; j<=knoten; j++) graph[i][j]=leer;
                            adj[i]=0;
                            drin[i][0]=drin[i][1]=0;
                            roots[i].tiefe=0;
                            blaetterkette[i]=nil;
                            blaetterzahl[i]=0; }
for (i=0; i<=knoten; i++) graph[0][i]=0;


}

/*****************************DEFAULTWERTE*****************************/

void defaultwerte(GRAPH graph, unsigned char *adj, unsigned char *startknoten)

//GRAPH graph;
//unsigned char *adj;
//unsigned char *startknoten;

/* schreibt die defaultwerte in einen zuvor initialisierten graphen und
   aktualisiert auch adj */

{

*startknoten=1;

graph[1][0]=2; graph[2][0]=1; stillness[1]-=2;
multiplicity=1;
adj[1]=adj[2]=1;
graph[0][0]=2;
graph[0][1]=1;
if (bipartite) { colour[1]=0; colour[2]=1; anzcolours[0]=anzcolours[1]=1;}


/* damit auch der minidar-baum vernuenftig initialisiert wird: */
if (treeborder>=1) minidar(graph,adj,(unsigned char)1);
}

/******************************GRAPHERZEUGUNG**************************/

void grapherzeugung()
{
GRAPH graph;
ADJAZENZ adj;
unsigned char startknoten;

initialize(graph,adj);

defaultwerte(graph,adj,&startknoten); 

v1_konstrukt(graph,adj,startknoten);

wegspeichern(globalliste,graphenzahl,graph,2); graphenzahl=0;
}

/*********************EINLESEN******************************************/


int einlesen(int argc, char *argv[])

//int argc;
//char *argv[];


{
int i, versch_valenzen=0;

if (argc < 2) 
    {
    fprintf(stderr,"Usage: multigraph n1 n2 ... [t<girth>] [b] [o|f] [m x y]\n");
    exit(1); 
    }

else { knotenzahl=kantenzahl=0; 
       for (i=0; i<=knoten; i++) orig_non[i]=nr_of_nodes[i]=0;
       for (i=1; (i<argc) && (argv[i][0] >= '0') && (argv[i][0]<='9') ; i++)
            { 
              orig_non[i]=nr_of_nodes[i]=atoi(argv[i]);
	      if (orig_non[i] != 0) versch_valenzen++;
              knotenzahl += nr_of_nodes[i];
              kantenzahl += i*nr_of_nodes[i]; }
       maxvalence=i-1;
       while ((maxvalence>0) &&(nr_of_nodes[maxvalence]==0)) maxvalence--;
       if (versch_valenzen==1) regular=1;

            if ((kantenzahl%2) || ((2*maxvalence)>kantenzahl) || 
                (knotenzahl>((kantenzahl/2)+1)) )
              {
                if (debug == 1) {
                    fprintf(stderr,"Solche zusammenhaengenden Graphen existieren nicht\n");
                }
                return 1;
              }
            else kantenzahl=kantenzahl/2;
            maxzykelkanten=kantenzahl-(knotenzahl-1); 
                         /*Die zyklomatische Zahl*/
            zykelkanten=0;
            taille=2;
       for ( ; i<argc; i++)
	   { switch(argv[i][0])
	     {
             case 'o': { output=1; noout=0; break; }
	     case 'f': { noout=0; break; }
	       /* in case both are used, standardout output dominates */
             case 'r': { recover=1; break; }
             case 't': { taille=atoi(argv[i]+1); break; }
             case 'B': { border=atoi(argv[i]+1); break; }
             case 'b': { bipartite=1; break; }
             case 'm': { i++; rest=atoi(argv[i]); i++; mod=atoi(argv[i]); break; }
             default: { fprintf(stderr,"%s nonidentified option\n",argv[i]); 
                        exit(0); }
	     }
	   }
       for (i=0; i<=knoten; i++) stillness[i]=0;
       for (i=maxvalence; i>=1; i--) 
                         stillness[i]=stillness[i+1]+nr_of_nodes[i];
       if (output && noout) { fprintf(stderr,"output on standardout AND just count?\n"); exit(1); }
     }
    return 0;
}

/**********************belegecomb****************************************/
void belegecomb(unsigned char a, unsigned char b, vertauschung *comb )

//unsigned char a;
//unsigned char b;
//vertauschung * comb;

/* belegt den vertauschungsvektor comb so, dass -- wenn diese vertauschungen
   auf einen vektor aus a+b Elementen angewendet werden, bei dem am Anfang
   elemente x0...x(a-1) stehen und am ende y0...y(b-1), alle permutationen
   erzeugt werden, die die reihenfolge von y0...y(b-1) respektieren -- nicht
   notwendig aber auch die von x0...x(a-1). Die Methode ist die von S.M.
   Johnson fuer Permutationen vorgeschlagene (Math.Comp.17 (1963) p.282-285,
   verallgemeinert wie von F.Ruskey vorgeschlagen (Congr.Num 67 (1988) 
   p.27-34). Die Implementierung ist vom Verfahren her umstaendlich und
   sicher nicht besonders effektiv -- das spielt aber keine Rolle, da der
   Laufzeitanteil am Programm verschwindend klein ist. */

{
  unsigned char positionen[knoten]; /* das array in dem die marken stehen */
  unsigned char I[knoten+1]; /* das I() aus Johnsons Arbeit */
  int i, zaehler;
  unsigned char wo, wohin, ende, kandidat;

  for(i=1; i<=a; i++) { positionen[i]=a-i+1; I[i]=0;}
  for(   ; i<= a+b; i++) { positionen[i]=0; I[i]=0; }

  comb[0][0]=comb[0][1]=0;
  zaehler=1;
  ende=0;
  while (!ende)
   { kandidat=0;
     for (i=1;i<=a+b;i++) 
        { if (positionen[i] > kandidat)
             { if (I[positionen[i]]==0)
                 { if ((i<a+b) && (positionen[i+1]<positionen[i]))
                      { kandidat=positionen[i]; wo=i; wohin=i+1; }
                 }
               else /* d.h. I[]=1 */
                 { if ((i>1) && (positionen[i-1]<positionen[i]))
                      { kandidat=positionen[i]; wo=i; wohin=i-1; }
                 }
             }
          } /*jetzt ist die naechste Vertauschung gefunden 
              -- oder kandidat==0 */
         if (kandidat==0) ende=1;
            else
            { for(i=kandidat+1; i<=a; i++) 
                        { if (I[i]) I[i]=0; else I[i]=1; }
              comb[zaehler][0]=wo; comb[zaehler][1]=wohin;
              zaehler++;
              positionen[wo]=positionen[wohin];
              positionen[wohin]=kandidat;
            }
     } /* ende while */
/*fprintf(stderr,"a:%d b: %d\n",a,b);
for(i=0;i<zaehler; i++) fprintf(stderr," %d %d \n",comb[i][0], comb[i][1]);
fprintf(stderr,"\n");*/
}


/*********************BELEGEVERT***************************************/

void belegevert(vertauschung vert[], int maxvalence)

//vertauschung vert[];
//int maxvalence;

{
    int i, j, k;

    vert[1][0] = 1;
    vert[1][1] = 2;

    for (i = 3; i <= maxvalence; i++)
    {
        for (j = 1; j < i; j++)
        {
            vert[fak[i - 1] + (j - 1) * fak[i - 1]][0] = i - j;
            vert[fak[i - 1] + (j - 1) * fak[i - 1]][1] = i - j + 1;
            /* der Wert n steht jetzt an der Stelle i-(j+1)*/

            for (k = 1; k <= fak[i - 1] - 1; k++)
            {
                vert[fak[i - 1] + (j - 1) * fak[i - 1] + k][0] = vert[k][0];
                vert[fak[i - 1] + (j - 1) * fak[i - 1] + k][1] = vert[k][1];
                if (vert[fak[i - 1] + (j - 1) * fak[i - 1] + k][0] >= i - j)
                    vert[fak[i - 1] + (j - 1) * fak[i - 1] + k][0]++;
                if (vert[fak[i - 1] + (j - 1) * fak[i - 1] + k][1] >= i - j)
                    vert[fak[i - 1] + (j - 1) * fak[i - 1] + k][1]++;
            }
        }
    }
}






/*********LAST*BUT*NOT*LEAST:****************MAIN***********************/

struct multigraph_result run_multigraph_gen(int argc, char** argv){
    int i, j;
    FILE* fil;
    char strpuf[60];
    struct tms TMS;


    recover = output = bipartite = regular = border = 0;
    /*Sicherungszeit=time(0)+Sicherungsintervall;*/

    int shouldExit = einlesen(argc, argv);
    if (shouldExit) {
        return (struct multigraph_result){NULL, -1, -1};
    }
    codelaenge = kantenzahl + knotenzahl;
    globalliste = (unsigned char*)malloc(codelaenge * listenlaenge);
    if (globalliste == NULL) {
        fprintf(stderr, "Can not get space for list of graphs \n");
        exit(0);
    }

    splitlevel = (kantenzahl * 3) / 5; /* nur geraten ....*/

    fak[0] = 1;
    fak[1] = 1;
    for (i = 2; i <= 8; i++) fak[i] = i * fak[i - 1];

    if (maxvalence > 8) {
        fprintf(stderr, "Die maximale Valenz ist zu gross !\n");
        exit(0);
    }

    for (i = 0; i < maxvalence; i++)
        for (j = 1; j <= maxvalence - i; j++) {
            comblaenge[i][j] = fak[i + j] / fak[j];
            combinations[i][j] =
                (vertauschung*)malloc(sizeof(vertauschung) * (comblaenge[i][j]));
            belegecomb(i, j, combinations[i][j]);
        }

    if (knotenzahl <= maxvalence) i = knotenzahl - 1;
    else i = maxvalence;
    vert = (vertauschung*)malloc(sizeof(vertauschung) * (fak[i]));
    if (vert == nil) {
        fprintf(stderr, "Nicht genug Platz zur Speicherung der Permutationen !\n");
        exit(0);
    }

    belegevert(vert, i);

    MASK[1] = 1;
    for (i = 2; i <= 32; i++) MASK[i] = MASK[i - 1] << 1;
    for (i = 33; i <= 64; i++) MASK[i] = MASK[i - 32];

    for (i = 1; i <= knoten; i++) {
        reihenfolge[i] = (unsigned char*)malloc(9 * sizeof(unsigned char));
        reihenfolge[i][0] = 0;
    }

    if (debug == 1) {
        for (i = 1; i <= maxvalence; i++)
            fprintf(stderr, "Knoten mit Grad %d : %d \n", i, nr_of_nodes[i]);
        fprintf(stderr, "maxvalence= %d; knotenzahl=%d\n", maxvalence, knotenzahl);
    }


    if (knotenzahl > knoten) {
        fprintf(stderr, "Knotenzahl zu gross !\n\n");
        exit(0);
    }

    graphenzahl = 0;

    sprintf(codename, "Codes");
    // This violates the C standard:
    //for (i=1;i<=maxvalence;i++) sprintf(codename,"%s.%d",codename,nr_of_nodes[i]);
    for (i = 1; i <= maxvalence; i++) sprintf(codename + strlen(codename), ".%d", nr_of_nodes[i]);

    if (taille >= 3) {
        for (i = argc - 1; argv[i][0] != 't'; i--);
        strcat(codename, argv[i]);
    }
    if (bipartite) strcat(codename, ".b");
    if (mod) {
        sprintf(strpuf, ".m_%d_%d", rest, mod);
        strcat(codename, strpuf);
    }
    if (!noout) {
        if (debug == 1) {
            if (!output) fprintf(stderr, "Output in das file %s !\n", codename);
            else fprintf(stderr, "Output auf stdout !\n");
        }
    }

    if ((!recover) && (!output) && !noout) {
        fil = fopen(codename, "w");
        fclose(fil);
    }

    /* Hier border sinnvoll ( wie auch immer ) festlegen */

    if (border == 0)
    {
        if (debug == 1) {
            fprintf(stderr, "BORDER muss noch festgelegt werden.\n");
        }
        border = kantenzahl - 8;
    }
    /*border=7;*/


    if (border - 2 <= baumgrenze) treeborder = border - 1;
    else treeborder = baumgrenze;

    if (treeborder < 0) treeborder = 0;

    if (debug == 1) {
        fprintf(stderr, "border=%d treeborder=%d \n", border, treeborder);
    }


    grapherzeugung();

    /*
          sprintf(strpuf,"rm ");
          strcat(strpuf,lastgraphname);

    system(strpuf);*/

    if (debug == 1) {
        fprintf(stderr, "Command:\n");
        for (i = 0; i < argc; i++) fprintf(stderr, "%s ", argv[i]);
        fprintf(stderr, "\n");
        times(&TMS);
#ifdef time_factor
        fprintf(stderr, "Erzeugung beendet. -- %.1f Sekunden.\n", (double)TMS.tms_utime / time_factor);
#else
        fprintf(stderr,"Erzeugung beendet. \n");
#endif
        fprintf(stderr, "Generated %d graphs with valence vector ", number_of_graphs);
        for (i = 1; i <= maxvalence; i++) fprintf(stderr, "%s ", argv[i]);
        fprintf(stderr, "and girth at least %d.\n", taille);
    }

    return (struct multigraph_result){globalliste, number_of_graphs, codelaenge};
}

int main_old(int argc, char* argv[]) {
    run_multigraph_gen(argc, argv);
    return 0;
}







