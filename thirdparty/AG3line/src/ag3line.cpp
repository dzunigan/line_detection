
#include "ag3line.h"

#include <ctime>
#include <queue>
#include <vector>
#include"BaseLine.h"
#include"GeneralFuncs.h"
#include"RegManager.h"
#include <cstdlib>
#include <string>
#include <fstream>
#include <ios>

/*----------------------------------------------------------------------------*/
/** Chained list of coordinates.
*/
struct coorlist
{
	ushort x, y;
	struct coorlist *next;
};
/*----------------------------------------------------------------------------*/
/*--------------------------- Rectangle structure ----------------------------*/
/*----------------------------------------------------------------------------*/


static void gaussian(const cv::Mat& im,cv::Mat& smooth)
{
	//gaussian filter
	//define the gaussian kernal
	float arr[]={0.00296901674395000,0.0133062098910137,0.0219382312797146,0.0133062098910137,0.00296901674395050,
		0.0133062098910137,0.0596342954361801,0.0983203313488458,0.0596342954361801,0.0133062098910137,
		0.0219382312797146,0.0983203313488458,0.162102821637127,0.0983203313488458,0.0219382312797146,
		0.0133062098910137,0.0596342954361801,0.0983203313488458,0.0596342954361801,0.0133062098910137,
		0.00296901674395050,0.0133062098910137,0.0219382312797146,0.0133062098910137,0.00296901674395050};
	cv::Mat kernel=cv::Mat(5,5,CV_32FC1,arr);
	filter2D( im, smooth, im.depth(),
		kernel);

}
static  void ang_grad(const cv::Mat& im,cv::Mat& grad,cv::Mat& ang,float gradt)
{
	//calculate the grad and angle data
	//define the x and y direction kernal
	int rsize=im.rows,colsize=im.cols;
	float dx[]={-1,1,-1,1};
	float dy[]={-1,-1,1,1};
	cv::Mat mdx=cv::Mat(2,2,CV_32FC1,dx);
	cv::Mat mdy=cv::Mat(2,2,CV_32FC1,dy);

	cv::Mat gx,gy;
	cv::filter2D(im, gx, im.depth(),mdx);
	cv::filter2D( im, gy, im.depth(),mdy);

	grad=cv::Mat::zeros(rsize, colsize, CV_32FC1);
	ang=cv::Mat::zeros(rsize, colsize, CV_32FC1);

	//access each pixel to calculate the angle
	float* ptrdx = (float*)gx.data;
	float* ptrdy = (float*)gy.data;
	float* ptrgrad = (float*)grad.data;
	float* ptrang = (float*)ang.data;
	float vdx,vdy,vgrad;
	for(int i = 0; i < rsize; i++)
	{
		for(int j = 0; j < colsize; j++)
		{
			vdx = *(ptrdx++);
			vdy = *(ptrdy++);
			vgrad=vdx*vdx+vdy*vdy;
			//reduce computation
			if (vgrad<GRADE_2)
			{
				*(ptrgrad++)=NOTDEF;
				*(ptrang++)=NOTDEF;

				continue;
			}

			*(ptrgrad++)=sqrt2(vgrad);

			//*(ptrang++)=cv::fastAtan2(vdy,vdx);
			*(ptrang++)=atan2approx(vdy,vdx);
			//cout<<std::atan2(vdy,vdx)<<" "<<atan2approx(vdy,vdx)<<endl;;
			//getchar();
		}
	}

}

static  struct coorlist* extractAnchor(cv::Mat ang,cv::Mat grad,cv::Mat& anchor)
{

	//extract the anchor
	int m=ang.rows;
	int n=ang.cols;
	int n_bins=1024;
	int idxbin,dirtype,id1,id2,id3,id4;
	float max_grad = 255.0;
	/* the rest of the variables are used for pseudo-ordering
	the gradient magnitude values */
	int list_count = 0;
	struct coorlist *list;
	struct coorlist **range_l_s;  /* array of pointers to start of bin list */
	struct coorlist **range_l_e;  /* array of pointers to end of bin list */
	struct coorlist *start;
	struct coorlist *end;


	/* get memory for "ordered" list of pixels */
	list = (struct coorlist *) std::calloc( (size_t) (m * n), sizeof(struct coorlist) );
	//*mem_p = (void *) list;
	range_l_s = (struct coorlist **) std::calloc( (size_t) n_bins,
		sizeof(struct coorlist *) );
	range_l_e = (struct coorlist **) std::calloc( (size_t) n_bins,
		sizeof(struct coorlist *) );

	for (int i = 0; i < n_bins; i++) range_l_s[i] = range_l_e[i] = NULL;



	anchor=cv::Mat::zeros(m,n,CV_8UC1);
	cv::Mat anchorrec=cv::Mat::zeros(m,n,CV_8UC1);
	int isanchor,isanchorgrad,binidx;
	float* ptrgrad = (float*)grad.data+n;
	float* ptrang = (float*)ang.data+n;
	uchar*ptranchor=(uchar*)anchor.data+n;
	uchar*ptrrec=(uchar*)anchorrec.data+n;
	//finish=clock();
	for(int i=1;i<m-1;i++)
	{
		for (int j=1;j<n-1;j++)
		{
			++ptrgrad;
			++ptrrec;
			++ptrang;
			++ptranchor;
			if (*ptrgrad==NOTDEF||*ptrrec==1)
				continue;

			dirtype=direction(*ptrang);



			switch (dirtype)
			{
			case 1:
				id1=1;
				id2=-1;
				id3=-n;
				id4=n;
				break;
			case 2:
				id1=n+1;
				id2=-n-1;
				id3=n-1;
				id4=-n+1;
				break;
			case 3:
				id1=-n;
				id2=n;
				id3=1;
				id4=-1;
				break;
			case 4:
				id1=n-1;
				id2=-n+1;
				id3=n+1;
				id4=-n-1;
				break;

			}




			if ((*ptrgrad)-*(ptrgrad+id1)<1.0||
				(*ptrgrad)-*(ptrgrad+id2)<1.0)
			{
				continue;
			}

			*(ptrrec+id1)=1;
			*(ptrrec+id2)=1;


			if (angle_diff(*ptrang,*(ptrang+id3))>PI_8||
				angle_diff(*ptrang,*(ptrang+id4))>PI_8)
			{
				*ptranchor=1;
				continue;
			}

			*ptranchor=2;

			/* store the Point3i in the right bin according to its norm */
			idxbin = (unsigned int) ((*ptrgrad) * (double) n_bins / max_grad);
			if ( idxbin >= n_bins ) idxbin = n_bins - 1;
			if ( range_l_e[idxbin] == NULL )
				range_l_s[idxbin] = range_l_e[idxbin] = list + list_count++;
			else
			{
				range_l_e[idxbin]->next = list + list_count;
				range_l_e[idxbin] = list + list_count++;
			}
			range_l_e[idxbin]->x = (ushort) j;
			range_l_e[idxbin]->y = (ushort) i;

			range_l_e[idxbin]->next = NULL;
		}
		//++ the address
		ptrgrad+=2;
		ptrrec+=2;
		ptrang+=2;
		ptranchor+=2;
	}


	/* Make the list of pixels (almost) ordered by norm value.
	It starts by the larger bin, so the list starts by the
	pixels with higher gradient value. Pixels would be ordered
	by norm value, up to a precision given by max_grad/n_bins.
	*/
	int i=0;
	for ( i= n_bins - 1; i > 0 && range_l_s[i] == NULL; i--);
	start = range_l_s[i];
	end = range_l_e[i];
	if ( start != NULL )
		for (i--; i > 0; i--)
		{
			if ( range_l_s[i] != NULL )
			{
				end->next = range_l_s[i];
				end = range_l_e[i];
			}
		}

	/* free memory */
	std::free( (void *) range_l_s );
	std::free( (void *) range_l_e );


	return start;
}

static void growLine(RegManager *rm,BaseLine *bl,uchar*anchor,float*ang,uchar*used,int maxitera)
{
	//Mat imc=imread("001.png",0);
	//imc=imc*0+255;
	int xx,yy,locat;
	bool aligned,found;
	int i,j;
	//Mat imc=Mat::zeros(ysize,xsize,CV_8UC1);
	float dist=0;
	while(1)
	{
		for(i=1;i<=maxitera;i++)
		{
			bl->newValidationPixel(i);
			found=false;
			// check if the three points are aligned
			for (j=0;j<3;j++)
			{
				xx=bl->xxs[j];
				yy=bl->yys[j];
				// check if the coor valid
				if(xx<0||xx>=bl->xsize||yy<0||yy>=bl->ysize)
					break;
				locat=yy*bl->xsize+xx;
				// check if this Point3i is anchor and have not been used
				if (used[locat]==USED||anchor[locat]==0)
					continue;

				// radius&anchor restrain
				if(i>3&&anchor[locat]==1)
					continue;
				//distance validation
				dist=bl->withinLength(xx,yy);
				if (dist>1.01)
					continue;


				//then the distance or angle restrain
				aligned=angle_diff(bl->segang,ang[locat])<=PI_8;
				if(!((bl->counter>bl->minlength&&i==1&&dist<=0.5)||aligned))
					continue;

				//add this Point3i
				rm->iteraarr[bl->counter]=i;
				rm->reg[bl->counter].x=xx;
				rm->reg[bl->counter].y=yy;
				rm->reg[bl->counter].z=locat;
				used[locat]=USED;

				//update the base line
				if (aligned)
				{
					bl->updateAnglePrams(ang[locat]);
				}

				//update the SVD parameters
				bl->updateLinePrams(xx,yy);
				found=true;
				break;


			}


			if (!found)
				continue;
			break;

		}
		if(!found)
			break;

	}

}

static void setUSED(cv::Point3i*reg,int testsize,uchar*used,int xsize,int ysize)
{
	int xx,yy;

	for(int i=0;i<testsize;i++)
	{
		xx=reg[i].x;
		yy=reg[i].y;
		if(xx<0||xx>=xsize||yy<0||yy>=ysize)
			continue;


		used[reg[i].z]=1;
		//cout<<"in"<<endl;
	}
}
static void getAnchors(cv::Point3i* rectreg,cv::Point3i*anchorreg,
	BaseLine*bl,int rectsize, int *anchorsize,uchar*anchor)
{

	*anchorsize=0;
	int xx,yy;
	for(int i=0;i<rectsize;i++)
	{
		xx=rectreg[i].x;
		yy=rectreg[i].y;

		//check if the location is right
		if(xx<0||xx>=bl->xsize||yy<0||yy>=bl->ysize)
			continue;
		if(anchor[rectreg[i].z]==0)
			continue;

		//aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa threshold
		//if(bl->withinLength(xx,yy)>1.414)
		//continue;
		*anchorsize=*anchorsize+1;
		anchorreg[*anchorsize].x=xx;
		anchorreg[*anchorsize].y=yy;
		anchorreg[*anchorsize].z=xx+yy*bl->xsize;


	}





}
static bool densityValidation(RegManager*rm,BaseLine*bl,uchar*used,uchar*anchor,float& thisdthre)
{
	int	linesize,anchorsize;

	int lastcounter=0,i;
	while (1)
	{
		//cannot split since the seg is unaccurate
		if (lastcounter==rm->counter)
			return false;
		lastcounter=rm->counter;
		//cout<<"gg1"<<endl;
		linesize=0;
		anchorsize=0;
		bl->region2line(rm->reg,rm->ed1,rm->ed2,rm->ed3,rm->ed4,PI_8);
		bl->getlength();

		//length validation
		if(!bl->lengthSatisfy())
		{
			//	cout<<"gg..............."<<endl;
			return false;
		}
		//linesize
		//calculate the line length

		EnumerateRectPoints1(bl->x1, bl->y1,bl->x2, bl->y2,rm->rectreg, &(rm->rectsize),2,bl->xsize);
		//anchorsize=getAnchorSize(rm->rectreg,testsize,bl,used,anchor);
		getAnchors(rm->rectreg,rm->anchorreg,bl,rm->rectsize,&anchorsize,anchor);

		thisdthre=(float)anchorsize/bl->length;

		//density validation
		if(thisdthre>=bl->getAnchorThreshold())
		{
			//setUSED(rm->anchorreg,anchorsize,anchor,used,bl->xsize,bl->ysize);
			return true;
		}

		//cannot satisfy the density
		//1st find the maxmum itera
		lastcounter=rm->counter;
		rm->splitTheReg();
	}

}
static void releaseUSED(uchar*used,cv::Point3i*reg,int ed1,int ed2,int ed3,int ed4,int m,int n,int mid,int total)
{
	int i;
	for( i=0;i<ed1;i++)
		used[reg[i].z]=0;
	for(i=ed2;i<mid;i++)
		used[reg[i].z]=0;
	for( i=mid;i<ed3;i++)
		used[reg[i].z]=0;
	for(i=ed4;i<total;i++)
		used[reg[i].z]=0;
}

static bool gradValidation(RegManager*rm,BaseLine*bl,uchar*used,uchar*anchor,float*grad)
{

	int linesize,anchorsize,idst,i,j,maxitera,maxidx,lastcounter=0;
	float std,alignment;

	//generate the validation grades
	while (1)
	{
		if (lastcounter==rm->counter)
			return false;

		EnumerateRectPoints1(bl->x1, bl->y1,bl->x2, bl->y2,rm->rectreg, &(rm->rectsize),1,bl->xsize);

		getAnchors(rm->rectreg,rm->anchorreg,bl,rm->rectsize,&anchorsize,anchor);

		//grade validation
		//1st get the grade from the linereg
		for( i=0;i<rm->rectsize;i++)
		{
			if (rm->rectreg[i].y>bl->ysize||rm->rectreg[i].x>bl->xsize
				||rm->rectreg[i].y<0||rm->rectreg[i].x<0)
				continue;

			rm->gradarr[i]=grad[rm->rectreg[i].z];
			if(rm->gradarr[i]==-1024)
				rm->gradarr[i]=5.2;
		}


		//calculate the parameter
		float m=0;
		sdeviation(rm->gradarr,rm->rectsize,&m,&std);
		alignment=m*(1-6*std/360);
		//plotGrad(rm,rm->rectsize,m,alignment);


		bool goodseg=true;
		int uncounter=0;
		//check every pixel on the line
		//cout<<rm->ed1<<" "<<rm->ed2<<" "<<rm->ed3<<" "<<rm->ed4<<endl;
		for (i=0;i<rm->rectsize;i++)
		{
			int xx=rm->rectreg[i].x;
			int yy=rm->rectreg[i].y;
			if(abs(grad[rm->rectreg[i].z]-m)<=alignment)
				continue;

			//check the density if there is a anchor within 3 distance(considering the fitting error) it can be viewd as a valid pixel
			//aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa threshold
			if (rm->satisfyGradAligned(anchorsize,i,3))
				continue;

			//add as a uncounter and check
			if (++uncounter>1)
			{
				goodseg=false;
				break;
			}
		}
		//return as a good line seg
		if (goodseg==true){
			return true;
		}
		lastcounter=rm->counter;
		rm->splitTheReg();
		if(rm->counter<bl->minlength)
			return false;
		//then refit the lineseg
		bl->region2line(rm->reg,rm->ed1,rm->ed2,rm->ed3,rm->ed4,PI_8);
		bl->getlength();

		//exit when the length is not satisfied
		if(!bl->lengthSatisfy())
			return false;
	}


}

bool firstThree(RegManager *rm,BaseLine *bl,uchar*anchor,float*ang,uchar*used,int maxitera)
{//check if the first three pixels satisfy the grouping

	int xx,yy;
	bool aligned,found;
	int i=2,j,locat;
	//Mat imc=Mat::zeros(ysize,xsize,CV_8UC1);
	float dist=0;
	while(i)
	{
		bl->newValidationPixel(1);
		found=false;
		// check if the three points are aligned
		for (j=0;j<3;j++)
		{
			xx=bl->xxs[j];
			yy=bl->yys[j];
			// check if the coor valid
			if(xx<0||xx>=bl->xsize||yy<0||yy>=bl->ysize)
				break;
			locat=yy*bl->xsize+xx;
			// check if this Point3i is anchor and have not been used
			if (used[locat]==USED||anchor[locat]<=1)
				continue;


			//distance validation
			dist=bl->withinLength(xx,yy);
			if (dist>0.6)
				continue;

			// angle restrain
			aligned=angle_diff(bl->segang,ang[locat])<=PI_16;
			if(!(aligned))
				continue;

			//add this Point3i
			rm->iteraarr[bl->counter]=1;
			rm->reg[bl->counter].x=xx;
			rm->reg[bl->counter].y=yy;
			rm->reg[bl->counter].z=locat;
			used[locat]=USED;

			//update the base line

			bl->updateAnglePrams(ang[locat]);


			//update the SVD parameters
			bl->updateLinePrams(xx,yy);
			found=true;
			break;
		}

		i--;
		if (!found)
			return false;

	}

	return true;
}

void ag3line(const cv::Mat& image,std::vector<lineag>& lines,bool control)
{
	cv::Mat im;
	image.convertTo(im, CV_32FC1);

	float gradt=5.2;//gradient threshold
	int initialSize=3;//the threshold of the initial short seg
	float ang_th = 22.5;     /* Gradient angle tolerance in degrees.           */
	float prec = M_PI * ang_th / 180.0;/* angle tolerance */
	float p = ang_th / 180.0;
	int maxitera=10;//the jump length
	int*binarr;

	int m,n;
	m=im.rows;
	n=im.cols;
	float N=m;
	if(n<m)
		N=n;

	int minlength=-4*std::log(N)/std::log(0.125);
	int minregsize=minlength*0.9;


	cv::Mat imsmoothM,gradM,angM,anchorM,usedM;
	struct coorlist *list_p;


	//surpres the noise
	gaussian(im, imsmoothM);

	//compute the gradient, angle and anchor
	ang_grad(imsmoothM,gradM,angM,gradt);

	list_p=extractAnchor(angM,gradM,anchorM);
	//drawAnchorMap(anchorM);

	usedM=cv::Mat::zeros(m,n,CV_8UC1);

	//initial some datas
	float* grad = (float*)gradM.data;
	float* ang = (float*)angM.data;
	uchar* used = (uchar*)usedM.data;
	uchar* anchor = (uchar*)anchorM.data;


	int x,y;
	float reg_angle,sumdx,sumdy,density;
	bool initialed=false;

	RegManager rm(m,n);
	BaseLine bl(m,n,minlength,117.8,-2.465,0.6996);
	bl.xsize=n;
	bl.ysize=m;
	int regsize,linesize;


	lineag lg;
	int cx,cy;
	//Mat im1=imread(picstr,0);
	//im1=im1*0;
	// ofstream pout(fname,ios::app);
	//loop the anchor for extraction
	for (; list_p != NULL; list_p = list_p->next )
	{
		if (!(used[ list_p->x + list_p->y * n ] == NOTUSED &&
			ang[list_p->x + list_p->y * n] != NOTDEF ))
			continue;//skip the used and cannot used pixels
		cx=list_p->x;
		cy=list_p->y;
		//the level 1 anchor cannot be used as a initial Point3i

		//initial some data
		rm.iteraarr[0]=1;
		rm.reg[0].x=cx;
		rm.reg[0].y=cy;
		rm.reg[0].z=cy*n+cx;
		bl.initial(cx,cy,ang[rm.reg[0].z]);

		///*
		if (!firstThree(&rm,&bl,anchor,ang,used,maxitera)){
			//release used
			for(int i=0;i<bl.counter;i++)
				used[rm.reg[i].z]=0;
			continue;
		}
		//*/
		//grow 2 the angle\s vertical direction
		growLine(&rm,&bl,anchor,ang,used,maxitera);
		rm.ed1=0;
		rm.ed2=bl.counter;
		rm.ed3=bl.counter;
		bl.midcounter=bl.counter;

		//reverse the direction and grow to the other side
		bl.reverseDir(cx,cy);
		growLine(&rm,&bl,anchor,ang,used,maxitera);
		rm.ed4=bl.counter;


		if(rm.ed4-rm.ed3+rm.ed2-rm.ed1<minlength)
			continue;
		if(!control)
		{
			bl.region2line(rm.reg,rm.ed1,rm.ed2,rm.ed3,rm.ed4,PI_8);
			EnumerateRectPoints1(bl.x1, bl.y1,bl.x2, bl.y2,rm.rectreg, &(rm.rectsize),2,bl.xsize);
		    setUSED(rm.rectreg,rm.rectsize,used,bl.xsize,bl.ysize);
			lg.x1=bl.x1;
			lg.y1=bl.y1;
			lg.x2=bl.x2;
			lg.y2=bl.y2;
			lines.push_back(lg);
			continue;
		}

		if (!densityValidation(&rm,&bl,used,anchor,density))
		{
			releaseUSED(used,rm.reg,rm.ed1,rm.ed2,rm.ed3,rm.ed4, m, n,bl.midcounter,bl.counter);
			continue;
		}


		// if the density is significant, save time
		if(density>0.98)
		{
			lg.x1=bl.x1;
			lg.y1=bl.y1;
			lg.x2=bl.x2;
			lg.y2=bl.y2;
			lines.push_back(lg);
			setUSED(rm.rectreg,rm.rectsize,used,bl.xsize,bl.ysize);

			continue;

		}
		// if the density is significant, save time

		if(!gradValidation(&rm,&bl,used,anchor,grad)){
			releaseUSED(used,rm.reg,rm.ed1,rm.ed2,rm.ed3,rm.ed4, m, n,bl.midcounter,bl.counter);
			continue;
		}
		//EnumerateRectPoints1(bl.x1, bl.y1,bl.x2, bl.y2,rm.rectreg, &(rm.rectsize),2,bl.xsize);

		setUSED(rm.rectreg,rm.rectsize,used,bl.xsize,bl.ysize);

		releaseUSED(used,rm.reg,rm.ed1,rm.ed2,rm.ed3,rm.ed4, m, n,bl.midcounter,bl.counter);


		lg.x1=bl.x1;
		lg.y1=bl.y1;
		lg.x2=bl.x2;
		lg.y2=bl.y2;
		lines.push_back(lg);

	}

	//pout.close();
	///*
	std::free( (void *) list_p );
}

/*
int main()
{

		bool control=true;//81550_150560154_042_20160205045406_rgb.jpg
		std::string imgstr="tower11.jpg";
		cv::Mat img = cv::imread(imgstr,0);
		cv::Mat imc = cv::imread(imgstr,0);
		cv::Mat im;




		std::vector<lineag> lines;
		std::clock_t start, finish;
		img.copyTo(im);
		start = std::clock();
		ag3line(im,lines,control);
		finish = std::clock();
		double duration = 1e3 * (finish - start) / CLOCKS_PER_SEC;
		//cout<<duration;
		//getchar();
		std::cout << duration << std::endl;
		imc=imc*0+255;
		for(int ii=0;ii<lines.size();ii++)
		{
			cv::line(imc,cv::Point(lines[ii].x1,lines[ii].y1),
				cv::Point(lines[ii].x2,lines[ii].y2),
				cv::Scalar(0,0,0),1,16);
		}
		cv::imwrite("resag3.jpg",imc);
		cv::imshow("ss",imc);
		cv::waitKey();

	return 0;
}
*/
