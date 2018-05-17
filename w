


#include <iostream>
#include <conio.h>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <windows.h>
#include <stdio.h>
#include <ctime>

using namespace std;



struct vector2
{ 
public: float x;
public: float y;

public: vector2(void) : x(0), y(0) 
{ 
}

public: vector2(float x1, float y1) : x(x1), y(y1) 
{ 
}

public: vector2(const vector2 & v) : x(v.x), y(v.y) 
{ 
}

public: vector2(const vector2 * v) : x(v->x), y(v->y) 
{ 
}


public: ~vector2(void) 
{ 
}


public: void set(float x1, float y1) 
{
x = x1;
y = y1;
}

public: vector2 maxi(vector2 a, vector2 b)
		{if (a.length() > b.length())
		{
		return a;
		}
		else 
		{
			return b;
		}
		}

public: float length() 
{
return sqrt(x * x + y * y);
}
		

public: float length_squared() 
{
return x * x + y * y;
}

public: float distance(const vector2 & v) 
{ 
return sqrt(((x - v.x) * (x - v.x)) + ((y - v.y) * (y - v.y)));
}

public: vector2 normal()
{ 
set(-y, x); 
return *this; 
}


public: vector2 normalize()
{

if(length() != 0)
{
    float length1 = sqrt(x*x + y*y);
    x = x/length1; 
	y = y/length1;
    return *this;
}
else 
{
	x = 0;
y = 0;
return *this;
}
}





public: vector2 operator = (const vector2 & v) 
{
x = v.x; 
y = v.y;
return *this;
}

public: vector2 operator = (const float & f) 
{
x = f;
y = f;
return *this;
}

public: vector2 operator - (void) 
{ 
x = -x; 
y = -y;
return *this;
}

public: bool vector2::operator == (const vector2 & v) const 
{ 
	return (x == v.x) && (y == v.y);
}

public: bool vector2::operator != (const vector2 & v) const 
{ 
	return (x != v.x) || (y != v.y);
}


public: vector2 operator + (const vector2 & v) const 
{
return vector2(x + v.x, y + v.y); 
}

public: vector2 operator - (const vector2 & v) const 
{ 
return vector2(x - v.x, y - v.y);
}



public: vector2 operator * (const vector2 & v) const 
{
return vector2(x * v.x, y * v.y); 
}

public: vector2 operator / (const vector2 & v) const 
{ 
return vector2(x / v.x, y / v.y);
}


public: vector2 operator += (const vector2 & v) 
{ 
x += v.x; 
y += v.y; 
return *this; 
}


public: vector2 operator -= (const vector2 & v) 
{ 
x -= v.x; 
y -= v.y;
return *this; 
}

public: vector2 operator *= (const vector2 & v) 
{ 
x *= v.x;
y *= v.y;
return *this; 
}

public: vector2  operator /= (const vector2 & v) 
{ 
x /= v.x;
y /= v.y; 
return *this;
}

public: vector2 operator + (float v) const 
{
return vector2(x + v, y + v);
}

public: vector2 operator - (float v) const
{ 
return vector2(x - v, y - v);
}

public: vector2 operator * (float v) const 
{
return vector2(x * v, y * v);
}

public: vector2 operator / (float v) const
{
return vector2(x / v, y / v);
}


public: vector2 operator += (float v) 
{ 
x += v; 
y += v; 
return *this;
}

public: vector2 operator -= (float v) 
{ 
x -= v;
y -= v; 
return *this;
}

public: vector2 operator *= (float v) 
{
x *= v;
y *= v; 
return *this;
}

public: vector2 operator /= (float v) 
{x /= v;
y /= v;
return *this;
}

};


class circle {
public:
	int centerx;
	int centery;
	int rad;
	bool nullbrush;
	HPEN pen;
	HBRUSH brush;


	void pencolor(COLORREF color)
	{if (pen)
	{
	DeleteObject(pen);
	}
	pen = CreatePen(PS_SOLID, 3, color);
	}


	void brushcolor(COLORREF color, bool Nullbrush = false)
	{
		nullbrush = Nullbrush;
		if (brush)
	{
		DeleteObject(brush);
	}

	if (Nullbrush)
	{ 
		brush = static_cast<HBRUSH>(GetStockObject(NULL_BRUSH));
		return;
	}

	brush = CreateSolidBrush(color);
	}


	
	circle (int x1, int y1, int r, bool Nullbrush = false) : centerx(x1), centery(y1), rad(r), pen(NULL), brush(NULL)
	{
	pencolor(RGB(0, 0, 0));
	brushcolor(RGB(0, 0, 0), Nullbrush);
	}

	circle () : centerx(0), centery(0), rad(0), pen(NULL), brush(NULL)
	{
	pencolor(RGB(0, 0, 0));
	brushcolor(RGB(0, 0, 0), false);
	}

	circle (int x1, int y1, int r, COLORREF penscolor, COLORREF paintscolor, bool Nullbrush = false) : centerx(x1), centery(y1), rad(r), pen(NULL), brush(NULL)
	{
	pencolor(penscolor);
	brushcolor(paintscolor, Nullbrush);
	}

	~circle()
	{
		if(pen)
		{
			DeleteObject(pen);
		}
		if (brush)
		{
			DeleteObject(brush);
	}
	}
	void draw(HDC hdc)
	{
		HPEN holdedpen = static_cast<HPEN>(SelectObject(hdc, pen));
        HBRUSH holdedbrush = static_cast<HBRUSH>(SelectObject(hdc, brush));
 
        Ellipse(hdc, centerx-rad, centery-rad, centerx+rad, centery+rad);
 
        SelectObject(hdc, holdedpen);
        SelectObject(hdc, holdedbrush);
	}

};

class colors
{public:
	
	int r;
int g;
int b;

colors (int r1, int g1, int b1) : r(r1),  g(g1), b(b1)
	{
	}

colors () : r(rand()%256),  g(rand()%256), b(rand()%256)
	{
	}

};

struct particle	
{
public: vector2 position;
public: vector2 dp;
public: float radius;
public: colors color;


		public: particle(void) : position(), dp(),  radius((float)0), color()
{ 
}

public: particle(float x1, float y1, float x2, float y2, float r, colors color1) : position(x1, x2), dp(x2, y2),  radius(r), color(color1)
{ 
		}

		public: particle(float x1, float y1, float x2, float y2, float r) : position(x1, x2), dp(x2, y2),  radius(r), color()
{ 
		}


        
	public: void Move()
        {
            position.x = dp.x + position.x;
			 position.y = dp.y + position.y;
            if (dp.length() > 3) 
            {
                dp.normalize();
                dp = dp*(float)3;
            }
        }

	public: void Push(vector2 delta)
        {
            position =  position + delta;
             dp = dp + delta;
        }
    };

   struct rectangle
   { public:
   float x;
   float y;
   float width;
   float height;



   public: rectangle(void) : width(0), height(0), x(0), y(0)
{ 
}
		   public: rectangle(float x2, float y2, float x1, float y1) : width(x2), height(y2), x(x1), y(y1)
{ 
}

   };


float randomfloat(float min, float max) {
    float random = ((float) rand()) / (float) RAND_MAX;
    float d = max - min;
    float r = random * d;
    return min + r;
}


struct link
{  public:
	float standartlenght;
	float stiffness;
	particle *particle1;
	particle *particle2;
	bool existance;

	link (float sl, float s, particle *q, particle *w, bool exist) 
	{
	standartlenght = sl;
	stiffness = s;
	*particle1 = *q;
	*particle2 = *w;
	existance = exist;
	}

	link () 
	{
	standartlenght = 0;
	stiffness = 0;
	particle *particle1 = new particle;
	particle *particle2 = new particle;
	existance = FALSE;
	}

	
link operator = (const link &A)
{
	this->standartlenght = A.standartlenght;
this->stiffness = A.stiffness;
this->particle1 = A.particle1;
this->particle2 = A.particle2;
this->existance = A.existance;
return *this;
}


	void pinkius_paius()
	{  
		if (stiffness == (float)1)
	{
		if ((particle1->dp.x)*(particle1->dp.x) > (particle2->dp.x)*(particle2->dp.x))
		{
			particle2->dp.x = particle1->dp.x;
		}
		else
		{
			particle1->dp.x = particle2->dp.x;
		}
	}

	else {

		vector2 p1 = particle1->position;
		vector2 p2 = particle2->position;
		vector2 norma = (p2 - p1).normalize();
		vector2 target1 = (p1 + p2) * (float)0.5 - norma *  standartlenght * (float)0.5;
		vector2 target2 = (p1 + p2) * (float)0.5 + norma *  standartlenght * (float)0.5;		
  particle1->position = particle1->position + (target1 - particle1->position) * stiffness;
  particle2->position = particle2->position + (target2 - particle2->position) * stiffness;
	}
	}


	void paint()
	{
		HDC hDC = GetDC( GetConsoleWindow( ) );
    HPEN Pen = CreatePen( PS_SOLID, 2, RGB(255, 255, 255));
    SelectObject( hDC, Pen );

	MoveToEx(hDC, particle1->position.x, particle1->position.y, NULL);
   LineTo(hDC, (int)particle2->position.x, (int)particle2->position.y);

	}
};




    class moreparticles
    {
	public: particle *particles;
        rectangle sandbox;
		link *links;


	public: moreparticles(int size, float radius, rectangle box, bool rnd)
        {
	    particles = new particle[size];
		links = new link[size+4];
	
            for (int i = 0; i < size ; i++)
            {
               
				if(rnd == TRUE)
                    {
						
					particles[i].position.x = randomfloat(box.x, box.width);
                    particles[i].position.y = randomfloat(box.y, box.height);
                    particles[i].radius = radius;
                    particles[i].dp = new vector2 (((float)rand()/(float)RAND_MAX - (float)0.5) * 100, ((float)rand()/(float)RAND_MAX - (float)0.5) * 100);
					
					
					bool a;
				    if (i == size+100)
					{
						a = TRUE;
					}
					else
					{
						a = FALSE;
					}
					if (i != size - 1)
					{
					links[i].standartlenght = radius*2;
					links[i].stiffness = (float)1;
					links[i].particle1 = (particles + i);
					links[i].particle2 = (particles + i+1); 
					links[i].existance = a;
					}
					if (i == size - 1)
					{
						links[i].existance = FALSE;
					}
				}



				if (rnd == false)
				{ 
				
					
				if(i < size/6 && 0 <= i)
					{
					particles[i].position.x = box.x + box.width/10 + 2*i*radius;
                    particles[i].position.y = box.y + box.height/4;
                    particles[i].radius = radius;
                    particles[i].dp = new vector2();

					if (i == 0)
					{
					links[i].standartlenght = (size/6-1)*radius*2;
					links[i].stiffness = (float)1;
					links[i].particle1 = (particles + i);
					links[i].particle2 = (particles + size/6 -1); 
					links[i].existance = TRUE;
					}
						if (i != 0)
					{
					links[i].standartlenght = radius*2;
					links[i].stiffness = (float)1;
					links[i].particle1 = (particles + i);
					links[i].particle2 = (particles + i - 1); 
					links[i].existance = TRUE;
					}
				}


				 if((i < 2*(size/6)) && (size/6 <= i))	
				 {
					particles[i].position.x = box.x + box.width/10 + 2*i*radius + 2*radius*(size/6);
                    particles[i].position.y = box.y + box.height/4;
                    particles[i].radius = radius;
                    particles[i].dp = new vector2();

					if (i == size/6)
					{
					links[i].standartlenght = (size/6-1)*radius*2;
					links[i].stiffness = (float)1;
					links[i].particle1 = (particles + i);
					links[i].particle2 = (particles + 2*(size/6) -1); 
					links[i].existance = TRUE;
					}
						if (i != size/6)
					{
					links[i].standartlenght = radius*2;
					links[i].stiffness = (float)1;
					links[i].particle1 = (particles + i);
					links[i].particle2 = (particles + i - 1); 
					links[i].existance = TRUE;
					}
				}

				 

				 if(i < 3*(size/6) && 2*(size/6) <= i)
					
				 {
					 particles[i].position.x = box.x + box.width/10 + (i - 2*(size/6))*radius*2;
                    particles[i].position.y = box.y + box.height/4 + 14*radius ;
                    particles[i].radius = radius;
					particles[i].dp = new vector2();

					if (i == (size/6)/5)
					{
					links[i].standartlenght = (size/6-1)*radius*2;
					links[i].stiffness = (float)1;
					links[i].particle1 = (particles + i);
					links[i].particle2 = (particles + 3*(size/6) -1); 
					links[i].existance = TRUE;
					}
						if (i != 2*(size/6))
					{
					links[i].standartlenght = radius*2;
					links[i].stiffness = (float)1;
					links[i].particle1 = (particles + i);
					links[i].particle2 = (particles + i - 1); 
					links[i].existance = TRUE;
					}
				}


					 if(i < 4*(size/6) && 3*(size/6) <= i)
					{
						particles[i].position.x = box.x + box.width/10 + (i - 2*(size/6))*radius*2+2*radius*(size/6);
                    particles[i].position.y = box.y + box.height/4+14*radius;
                    particles[i].radius = radius;
					particles[i].dp = new vector2();

					if (i == 3*(size/6))
					{
					links[i].standartlenght = (size/6-1)*radius*2;
					links[i].stiffness = (float)1;
					links[i].particle1 = (particles + i);
					links[i].particle2 = (particles + 4*(size/6) -1); 
					links[i].existance = TRUE;
					}
						if (i != 3*(size/6))
					{
					links[i].standartlenght = radius*2;
					links[i].stiffness = (float)1;
					links[i].particle1 = (particles + i);
					links[i].particle2 = (particles + i - 1); 
					links[i].existance = TRUE;
					}
					}

					 if(i < 6*(size/6) && 4*(size/6) <= i)
					{
						particles[i].position.x =  (particles[(size/6)*2+(size/6)/2].position.x + particles[(size/6)*2-(size/6)/2].position.x)/2 +  (i - 4*(size/6))*radius*3 - (size/6)*3*radius + 0.5*radius;
                    particles[i].position.y = box.y + box.height/4 + 7*radius;
                    particles[i].radius = radius;
					particles[i].dp = new vector2();
					

					if (i == 4*(size/6))
					{
					links[i].standartlenght = (size/6-1)*radius*3 + 0.5*radius;
					links[i].stiffness = (float)0;
					links[i].particle1 = (particles + i);
					links[i].particle2 = (particles + 5*(size/6) -1); 
					links[i].existance = TRUE;
					}
						if (i != 4*(size/6))
					{
					links[i].standartlenght = radius*3;
					links[i].stiffness = (float)0;
					links[i].particle1 = (particles + i);
					links[i].particle2 = (particles + i - 1); 
					links[i].existance = TRUE;
					}

						if (i == (size/6)*4)
					{
						particles[i].position.x = particles[i].position.x - 0.5*radius;
						
					}
						if (i == (size/6)*4+1)
						{
							links[i].standartlenght = radius*3.5;
						}

						if (i == (size/6)*6 - 1)
						{
							particles[i].position.x = particles[i].position.x + 0.5*radius;
							links[i].standartlenght = radius*3.5;
						}

					}



					 if(i == 6*(size/6) )
					 {
						 particles[i].position.x = particles[4*(size/6)].position.x ;
                    particles[i].position.y = particles[4*(size/6)].position.y + 2.5*radius;
                    particles[i].radius = radius;
					particles[i].dp = new vector2();
					links[i].standartlenght = radius*2.5;
					links[i].stiffness = (float)0;
					links[i].particle1 = (particles + i);
					links[i].particle2 = (particles + 4*(size/6)); 
					links[i].existance = TRUE;
					 }

					 if(i == 6*(size/6) + 1 )
					 {
						 particles[i].position.x = particles[6*(size/6)-1].position.x ;
                    particles[i].position.y = particles[6*(size/6)-1].position.y + 2.5*radius;
                    particles[i].radius = radius;
					particles[i].dp = new vector2();
					links[i].standartlenght = radius*2.5;
					links[i].stiffness = (float)0;
					links[i].particle1 = (particles + i);
					links[i].particle2 = (particles + 6*(size/6)-1); 
					links[i].existance = TRUE;
					 }

					  if(i == 6*(size/6) +2 )
					 {
						 particles[i].position.x = particles[4*(size/6)].position.x ;
                    particles[i].position.y = particles[4*(size/6)].position.y - 2.5*radius;
                    particles[i].radius = radius;
					particles[i].dp = new vector2();
					links[i].standartlenght = radius*2.5;
					links[i].stiffness = (float)0;
					links[i].particle1 = (particles + i);
					links[i].particle2 = (particles + 4*(size/6)); 
					links[i].existance = TRUE;
					 }

					 if(i == 6*(size/6)+ 3 )
					 {
						 particles[i].position.x = particles[6*(size/6)-1].position.x ;
                    particles[i].position.y = particles[6*(size/6)-1].position.y - 2.5*radius;
                    particles[i].radius = radius;
					particles[i].dp = new vector2();
					links[i].standartlenght = radius*2.5;
					links[i].stiffness = (float)0;
					links[i].particle1 = (particles + i);
					links[i].particle2 = (particles + 6*(size/6)-1); 
					links[i].existance = TRUE;
					 }
					
				}
										
                
            }

			
            sandbox = box;
        }

	public: void process(int size)
        {
			for (int i = 0; i < size + 4; i++)
            {   for (int j = 0; j < size+4; j++)
			{
				if (links[i].existance == TRUE)
				{  
						links[i].pinkius_paius();
				
					//links[i].paint();
					
				}
			}
            }
			for (int i = 0; i < size + 4; i++)
            {   for (int j = 0; j < size+4; j++)
			{
				if (links[i].existance == TRUE)
				{  
						links[i].pinkius_paius();
				
					//links[i].paint();
					
				}
			}
            }


            for (int i = 0; i < size - 1; i++)
             
				{
					for (int j = i + 1; j < size; j++)
                {
                    vector2 p1 = particles[i].position;
                    vector2 p2 = particles[j].position;
					
                    if ((p1 - p2).length_squared() < (particles[i].radius + particles[j].radius) * (particles[i].radius + particles[j].radius))
                    {

                        vector2 direction = p2 - p1;
                        direction.normalize();

                        float depth = particles[i].radius + particles[j].radius - (p2 - p1).length();

                        particles[i].Push(-(direction * depth * (float)0.5));
                       particles[j].Push(direction * depth * (float)0.5);

                      

                        vector2 relativevelocity = particles[i].dp - particles[j].dp;
                        
                        float velocity = ((float)1.5) * (relativevelocity.length() * direction.length());
                        if (velocity > 0)
                        {
                            particles[i].dp = particles[i].dp - direction * velocity *  (float)0.5;
                           particles[j].dp = particles[i].dp + direction * velocity *  (float)0.5;
						}
					

                    }
                }
			}

  
            for (int i = 0; i < size ; i++)
            {
                if (particles[i].position.x < sandbox.x + particles[i].radius)
                {
                    particles[i].position.x = sandbox.x + particles[i].radius;
                    if (particles[i].dp.x < (float)0.0)
                        particles[i].dp.x = (float)abs(particles[i].dp.x);// 

                }
                if (particles[i].position.x > sandbox.x + sandbox.width - particles[i].radius)
                {
                    particles[i].position.x = sandbox.x + sandbox.width - particles[i].radius;
                    if (particles[i].dp.x > (float)0.0)
                        particles[i].dp.x = -(float)abs(particles[i].dp.x);

                }
                if (particles[i].position.y < sandbox.y + particles[i].radius)
                {
                    particles[i].position.y = sandbox.y + particles[i].radius;
                    if (particles[i].dp.y < (float)0.0)
                        particles[i].dp.y = (float)abs(particles[i].dp.y);

                   
                }
                if (particles[i].position.y > sandbox.y + sandbox.height - particles[i].radius)
                {
                    particles[i].position.y = sandbox.y + sandbox.height - particles[i].radius;
                    if (particles[i].dp.y > (float)0.0)
                        particles[i].dp.y = -(float)abs(particles[i].dp.y);// 


                }
            }



			 
            for (int i = 0; i < size + 4; i++)
            {   
                particles[i].Move();
				
            }
        }

			public: void start(int size, rectangle box)
					{
						
						
						particles[6*(size/6)].dp.y = 1;
						particles[6*(size/6) + 1].dp.y = 1;
						particles[6*(size/6) + 2].dp.y = -1;
						particles[6*(size/6) + 3].dp.y = -1;
						
					}




					void process1(int size)
					{  
						 for (int i = 0; i < (size/6)*4 ; i++)
						 {
							 for (int j = 0; j < 4; j++)
							 { 
					vector2 p1 = particles[i].position;
                    vector2 p2 = particles[6*(size/6)+j].position;
					if ((p1 - p2).length_squared() <= ( (particles[i].radius + particles[6*(size/6)+j].radius) * (particles[i].radius + particles[6*(size/6)+j].radius)*2))
					{
						particles[6*(size/6) + j].dp.y = 0;
						links[size+j].standartlenght = particles[i].radius + particles[6*(size/6)+j].radius;
					links[size+j].stiffness = (float)1;
					links[size+j].particle1 = (particles + i);
					links[size+j].particle2 = (particles + 6*(size/6) + j); 
					links[size+j].existance = TRUE;
					links[6*(size/6)+j].standartlenght = (links[6*(size/6)+j].particle1->position-links[6*(size/6)+j].particle1->position).length();
					links[6*(size/6)+j].stiffness = (float)1;
					
					
					}

							 }
							 

						 
						 }
					
					}
					

  void process2(int size, rectangle box, int i)
           {  
        



			if( i < 5*(size/6) && 4*(size/6) <= i)
			{ 
				particles[i].dp.x = 0;
				particles[5*(size/6) + 5*(size/6)-i-1].dp.x = 0;
				if (-particles[i].position.x + particles[5*(size/6)+5*(size/6)-i-1].position.x > 2*(5*(size/6)-i)*(particles[i].radius +  particles[5*(size/6)+5*(size/6)-i-1].radius) - (particles[i].radius +  particles[5*(size/6)+5*(size/6)-i-1].radius) + 2*(5*(size/6)-i) )
			{
				
			 particles[i].dp.x =  0.05;
			 particles[5*(size/6) + 5*(size/6)-i-1].dp.x =  -0.05;
			 
			 
			 
				}
			
			 
  
  }
     



  }
			






public:	void paint(int size, circle *circles)
				{ 
					HWND hWnd = GetConsoleWindow();
    HDC hDC = GetDC(hWnd);
    RECT bmpRect;
    GetClientRect(hWnd,&bmpRect);
    HDC hBufferDC = CreateCompatibleDC(hDC);
    HBITMAP hBufferBmp=CreateBitmap(bmpRect.right,bmpRect.bottom,1,32,NULL);
    HBITMAP hBufferBmpOld=(HBITMAP)SelectObject(hBufferDC, hBufferBmp);


					for (int i = 0; i < size; i++)
				{
					
    circles[i].centerx = (int) particles[i].position.x;
	circles[i].centery = (int) particles[i].position.y;
	circles[i].rad = (int)particles[i].radius;
	circles[i].pencolor(RGB(particles[i].color.r, particles[i].color.g, particles[i].color.b));
	circles[i].brushcolor(RGB(particles[i].color.r, particles[i].color.g, particles[i].color.b));
	circles[i].draw(hBufferDC);
	BitBlt(hDC,0,0,bmpRect.right,bmpRect.bottom,hBufferDC,0,0,SRCCOPY);

				}
				//Sleep(5);
				}
				
    };



	int main()
	{   
			 srand(time(0));

  int z;
z = 11;

		 HWND hWnd = GetConsoleWindow();
    HDC hDC = GetDC(hWnd);
    RECT bmpRect;
    GetClientRect(hWnd,&bmpRect);
    HDC hBufferDC = CreateCompatibleDC(hDC);
    HBITMAP hBufferBmp=CreateBitmap(bmpRect.right,bmpRect.bottom,1,32,NULL);
    HBITMAP hBufferBmpOld=(HBITMAP)SelectObject(hBufferDC, hBufferBmp);
	

	
	wchar_t szTITLE[] = L"Pleasekillme";

	SetConsoleTitle(szTITLE);
	Sleep(1.5);
	MoveWindow(FindWindow(NULL, szTITLE), 0, 0, 600, 600, false);



 
  int a;
  a = 6*z+4;
  bool b;
  b = false;
  int k;
  k = 8;
		rectangle A((float)600, (float)600, (float)0, (float)0);
		moreparticles X(a, k, A, b);
		circle *circles = new circle[a];
		
		if (b == false)
		{ 
			X.start(a, A);
		}

		
    while(!GetAsyncKeyState(VK_RETURN))
		{
			X.paint(a, circles);
		
			

	   while (X.links[a].existance == FALSE && X.links[a+1].existance == FALSE && X.links[a+2].existance == FALSE && X.links[a+3].existance == FALSE && !GetAsyncKeyState(VK_RETURN) )
	   {
		  X.paint(a, circles);
		  
		X.process(a);
		X.process1(a);
		for (int o = 0; o < a+4; o++)
		{
			if(X.links[o].existance == TRUE)
			{
				X.links[o].paint();
			}
		}
		//Sleep(50);
	   }

	   if (b == false )
	   {  
		   while (X.particles[(a/6)].position.x - X.particles[(a/6)-1].position.x > 4*X.particles[(a/6)-1].radius)
	   {
	      X.paint(a, circles); 
		   for (int o = 0; o < a+4; o++)
		{
			if(X.links[o].existance == TRUE)
			{
				X.links[o].paint();
			}
		}

		   for (int h = 0; h < a + 4; h++)
		   {
			   X.process2(a, A, h);
			 
			
			
		 if ((X.particles[(a/6)].position.x - X.particles[(a/6)-1].position.x) > 4*X.particles[(a/6)-1].radius && - X.particles[(a/6)*4].position.x + X.particles[5*(a/6)+5*(a/6)-(a/6)*4-1].position.x > 2*(5*(a/6)-(a/6)*4)*(X.particles[(a/6)*4].radius +  X.particles[5*(a/6)+5*(a/6)-(a/6)*4-1].radius) - (X.particles[(a/6)*4].radius +  X.particles[5*(a/6)+5*(a/6)-(a/6)*4-1].radius) + 2*(5*(a/6)-(a/6)*4) )
		X.process(a);
			 
		
		
		   }
	   }
	   }
	   
	

		}

		Sleep (100000);
		
		return 0;
	}
