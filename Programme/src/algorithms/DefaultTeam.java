package algorithms;

import java.awt.Point;
import java.util.ArrayList;

import supportGUI.Circle;
import supportGUI.Line;

public class DefaultTeam {

  // calculDiametre: ArrayList<Point> --> Line
  //   renvoie une pair de points de la liste, de distance maximum.
  public Line calculDiametre(ArrayList<Point> points) {
    if (points.size()<3) {
      return null;
    }

    Point p=points.get(0);
    Point q=points.get(1);

    /*******************
     * PARTIE A ECRIRE *
     *******************/
    return new Line(p,q);
  }

  // calculDiametreOptimise: ArrayList<Point> --> Line
  //   renvoie une pair de points de la liste, de distance maximum.
  public Line calculDiametreOptimise(ArrayList<Point> points) {
    if (points.size()<3) {
      return null;
    }

    Point p=points.get(1);
    Point q=points.get(2);

    /*******************
     * PARTIE A ECRIRE *
     *******************/
    return new Line(p,q);
  }

  // calculCercleMin: ArrayList<Point> --> Circle
  //   renvoie un cercle couvrant tout point de la liste, de rayon minimum.
  public Circle calculCercleMin(ArrayList<Point> points) {
    if (points.isEmpty()) {
      return null;
    }

    Point center=points.get(0);
    int radius=100;

    /*******************
     * PARTIE A ECRIRE *
     *******************/
    return new Circle(center,radius);
  }

  // enveloppeConvexe: ArrayList<Point> --> ArrayList<Point>
  //   renvoie l'enveloppe convexe de la liste.
  
  public ArrayList<Point> enveloppeConvexe(ArrayList<Point> points){
    if (points.size()<3) {
      return null;
    }

    /*Algo naif:
    ArrayList<Point> enveloppe = new ArrayList<Point>();
     for(Point p : points) {
    	for(Point q : points) {
    		if(p.equals(q)) continue;
    		boolean estCote = true;
    		Point ref = p;
    		for(Point r : points) {
    			  if(crossProduct(p,q,p,r) != 0) {
    				  ref = r;
    				  break;
    			  }
    		  }
    		if(ref.equals(p)) {enveloppe.add(p);enveloppe.add(q);continue;}
    		for (Point r : points) {
    			if (crossProduct(p,q,p,ref)*crossProduct(p,q,p,r) < 0) {
    				estCote = false;
    				break;
    			}
    		}
    		if(estCote) {enveloppe.add(p);enveloppe.add(q);}
    	}
    }
    return enveloppe;*/
	
    Point A = new Point(0,0);
	Point B = new Point(Integer.MAX_VALUE,0);
	Point C = new Point(Integer.MAX_VALUE,Integer.MAX_VALUE);
	Point D = new Point(0,Integer.MAX_VALUE);
    
	//BIN SORT
    //on défini un tableau d'ymin avec une limite de 4000point pour la vitesse de traitement
    Point[] ymin = new Point[4000];
    for(Point p : points) {
    	if(ymin[p.x] == null || p.y < ymin[p.x].y) {
    		ymin[p.x]=p;
    	}
    }
    ArrayList<Point> parcResult = new ArrayList<Point>();
    for(int i=0; i < ymin.length; i++) {
    	if(ymin[i] != null) {
    		parcResult.add(ymin[i]);
    	}
    }
    
    Point[] ymax = new Point[4000];
    for(Point p : points) {
    	if(ymax[p.x] == null || p.y > ymax[p.x].y) {
    		ymax[p.x] = p;
    	}
    }
    for(int i = ymin.length-1 ; i >= 0 ; i--) {
    	if(ymax[i] != null) {
    		parcResult.add(ymax[i]);
    	}
    }
    //Algo de Graham
    double size = parcResult.size();
    for(int i=0; i < size + 16; i++) {
	   Point p = parcResult.get((int) (i%size));
	   Point q = parcResult.get((int) ((i+1)%size));
	   Point r = parcResult.get((int) ((i+2)%size));
	   if(crossProduct(p,q,p,r) < 0) {
		   parcResult.remove((int)((i+1)%size));
		   size = parcResult.size();
		   i = (int)((i-16)%size);
		   if(i<0) i=0;
	   }
    }
    //Algo naif pour recMin
    
    //test version
    
    
    ArrayList<Point> recMin = new ArrayList<Point>();
    for(int i=0; i < parcResult.size(); i++) {
    	Point p = parcResult.get(i);
    	Point q = parcResult.get((i+1)%parcResult.size());
    	Point s = p;
    	for(Point r: points) {
    		if(crossProduct(p,q,p,r) > crossProduct(p,q,p,s)) s=r;
    	}
    	double alphaPQ = (p.y-q.y)/(double)(p.x-q.x);
    	double betaPQ = p.y-alphaPQ*p.x;
    	double alphassPrime = 1/(double)alphaPQ;
    	double betassPrime = s.y-alphassPrime*s.x;
    	double sPrimeX = (betaPQ-betassPrime)/(alphassPrime-alphaPQ);
    	double sPrimeY = alphassPrime*sPrimeX+betassPrime;
    	Point sp = new Point(0,0);
    	sp.setLocation(sPrimeX, sPrimeY);
    	Point t = sp;
    	Point u = sp;
    	for(Point r: parcResult) {
    		if(crossProduct(sp,s,sp,r) > crossProduct(sp,s,sp,t)) t = r;
    		if(crossProduct(sp,s,sp,r) < crossProduct(sp,s,sp,u)) u = r;
    	}
    	double alphaGamma = alphassPrime;
    	double betaGamma = t.y-alphaGamma*t.x;
    	double alphaLambda = alphassPrime;
    	double betaLambda = u.y-alphaLambda*u.x;
    	double alphaPhi = alphaPQ;
    	double betaPhi = s.y-alphaPhi*s.x;
    	double aPrimeX = (betaPQ-betaLambda)/(alphaLambda-alphaPQ);
    	double aPrimeY = alphaPQ*aPrimeX+betaPQ;
    	double bPrimeX = (betaPQ-betaGamma)/(alphaGamma-alphaPQ);
    	double bPrimeY = alphaPQ*bPrimeX+betaPQ;
    	double cPrimeX = (betaPhi-betaGamma)/(alphaGamma-alphaPhi);
    	double cPrimeY = alphaPhi*cPrimeX+betaPhi;
    	double dPrimeX = (betaLambda-betaPhi)/(alphaPhi-alphaLambda);
    	double dPrimeY = alphaLambda*dPrimeX+betaLambda;
    	double aire = Math.sqrt(Math.pow(B.x-A.x,2)+Math.pow(B.y-A.y,2))*Math.sqrt(Math.pow(C.x-A.x,2)+Math.pow(C.y-A.y,2));
    	double airePrime = Math.sqrt(Math.pow(bPrimeX-aPrimeX,2)+Math.pow(bPrimeY-aPrimeY,2))*Math.sqrt(Math.pow(cPrimeX-aPrimeX,2)*Math.pow(cPrimeY-aPrimeY,2));
    	
    	if (airePrime < aire) {
    		A.setLocation(aPrimeX, aPrimeY);
    		B.setLocation(bPrimeX, bPrimeY);
    		C.setLocation(cPrimeX,cPrimeY);
    		D.setLocation(dPrimeX,dPrimeY);
    		recMin.clear();
    		recMin.add(A);
    		recMin.add(B);
    		recMin.add(C);
    		recMin.add(D);
    	}
    	
    }
    return recMin;
  } 
    
    /*for (int i = 0; i=< (enveloppe.size(); i++){
        Point A = parcResult.get (i);
        Point B = parcResult.get ((i+1)% parcResult.size());
    //1.3
        Point s = p;
        for (Point r: parcResult) if (crossProduct(p,q,p,r) > crossProduct(p,q,p,s)) s = r;
    //1.4 droite (s,s')
        double alphapq = (p.y - q.y)/(double)(p.x - q.x);
        double betapq = p.y - alphapq * p.x;
        double alphassPrime = 1/(double)alphapq;
        double betassPrime = s.y - alphassPrime * s.x;
        double sPrimeX = (betapq - betassPrime)/(alphassPrime - alphapq);
        double sPrimeY = alphassPrime * sPrimeX + betassPrime;
    //1.5
        Point sP = new Point(0,0);
        sP.setlocation(sPrimeX, sPrimeY);
        Point t = sP;
        Point u = sP;
        for (Point r: parcResult){
            if (crossProduct(sP,s,sP,r) > crossProduct(sP,s,sP,t)) t = r;
            if (crossProduct(sP,s,sP,r) > crossProduct(sP,s,sP,u)) t = u;
        }
    //1.6
        double alphaGamma = alphassPrime;
        double betaGamma = t.y - alphaGamma * t.x ;
    //1.7
        double alphaLambda = alphassPrime;
        double betaLambda = u.y - alphaLambda * u.x ;
    //1.8
        double alphaPhi = alphapq;
        double betaPhi = s.y - alphaPhi * s.x;
    //1.9
        Point a = new Point(0,0);
        Point b = new Point(0,Integer.Max.int);
        Point c = new Point(Integer.Max.int,Integer.Max.int);
        Point d = new Point(Integer.Max.int,0);

        double aPrimeX = (betapq - betaLambda)/(alphaLambda - alphapq);
        double aPrimeY = (alphapq * aPrimeX + betapq);

        double bPrimeX = (betapq - betaGamma)/(alphaGamma - alphapq);
        double bPrimeY = (alphapq * bPrimeX + betapq);

        double cPrimeX = (betaPhi - betaGamma)/(alphaGamma - alphaPhi);
        double cPrimeY = (alphaPhi * cPrimeX + betaPhi);

        double dPrimeX = (betaLambda - betaPhi)/(alphaPhi - alphaLambda);
        double dPrimeY = (alphaLambda * dPrimeX + betaLambda);
    //1.10
        double aire = Math.sqrt(Math.pow(B.x-A.x,2)+Math.pow(B.y-A.y,2))*Math.sqrt(Math.pow(C.x-A.x,2)+Math.pow(C.y-A.y,2));
        double airePrime = Math.sqrt(Math.pow(bPrimeX - aPrimeX,2) + Math.pow(bPrimeY - aPrimeY,2)) * Math.sqrt(Math.pow(cPrimeX - aPrimeX,2) + Math.pow(cPrimeY - aPrimeY,2));
            if (airePrime < aire) {
                A.setLocation(aPrimeX, aPrimeY);
                B.setLocation(bPrimeX, bPrimeY);
                C.setLocation(cPrimeX,cPrimeY);
                D.setLocation(dPrimeX,dPrimeY);
                recMin.clear();
                recMin.add(A);
                recMin.add(B);
                recMin.add(C);
                recMin.add(D);
        	}
        	
        }
        return recMin;
    }*/
    
    
  //Produit vectoriel
  private double crossProduct(Point p, Point q, Point s, Point t) {
	  return (q.x-p.x)*(t.y-s.y)-(q.y-p.y)*(t.x-s.x);
  }
 
    
  
}
