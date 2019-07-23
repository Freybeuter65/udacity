/* \author Aaron Brown */
// Quiz on implementing simple RANSAC line fitting

#include <cmath>
#include "../../render/render.h"
#include <unordered_set>
#include "../../processPointClouds.h"
// using templates for processPointClouds so also include .cpp to help linker
#include "../../processPointClouds.cpp"


pcl::PointCloud<pcl::PointXYZ>::Ptr CreateData()
{
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>());
  	// Add inliers
  	float scatter = 0.6;
  	for(int i = -5; i < 5; i++)
  	{
  		double rx = 2*(((double) rand() / (RAND_MAX))-0.5);
  		double ry = 2*(((double) rand() / (RAND_MAX))-0.5);
  		pcl::PointXYZ point;
  		point.x = i+scatter*rx;
  		point.y = i+scatter*ry;
  		point.z = 0;

  		cloud->points.push_back(point);
  	}
  	// Add outliers
  	int numOutliers = 10;
  	while(numOutliers--)
  	{
  		double rx = 2*(((double) rand() / (RAND_MAX))-0.5);
  		double ry = 2*(((double) rand() / (RAND_MAX))-0.5);
  		pcl::PointXYZ point;
  		point.x = 5*rx;
  		point.y = 5*ry;
  		point.z = 0;

  		cloud->points.push_back(point);

  	}
  	cloud->width = cloud->points.size();
  	cloud->height = 1;

  	return cloud;

}


pcl::PointCloud<pcl::PointXYZ>::Ptr CreateData3D()
{
	ProcessPointClouds<pcl::PointXYZ> pointProcessor;
	return pointProcessor.loadPcd("../../../sensors/data/pcd/simpleHighway.pcd");
}


pcl::visualization::PCLVisualizer::Ptr initScene()
{
	pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer ("2D Viewer"));
	viewer->setBackgroundColor (0, 0, 0);
  	viewer->initCameraParameters();
  	viewer->setCameraPosition(0, 0, 15, 0, 1, 0);
  	viewer->addCoordinateSystem (1.0);
  	return viewer;
}


std::unordered_set<int> RansacPlane(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, int maxIterations, float distanceTol )
{
	std::unordered_set<int> inliersResult[2];
	srand(time(NULL));

	double dist, dist_cnt, dist_sqrt;
	struct
	{
		float pln_cff_a;
		float pln_cff_b;
		float pln_cff_c;
		float pln_cff_d;
	} qtn;	///equation with the four coefficients
	int p1, p2, p3;
	int inlier_idx = 0;
	int inlier_out = 0;


	std::cout << "Cloud_Size: " << "\t" << cloud->points.size() << "\n";

	int iter = maxIterations;
	while( --iter >= 0  )
	{
	 	p1 = rand() % 976;
	  	p2 = rand() % 976;
	  	p3 = rand() % 976;


	// for( p1 = 0; p1 < cloud->points.size(); p1++ )
	// {
	// 	int p2 = p1 + 1;
	//  	for( ; p2 < cloud->points.size(); p2++ )
	//  	{
	//  		int p3 = p2 + 1;
	//  		for( ; p3 < cloud->points.size(); p3++ )
	//  		{
				// get the plane points pnt1, pnt2, pnt3
				pcl::PointXYZ pnt1 = cloud->points[p1];
				pcl::PointXYZ pnt2 = cloud->points[p2];
				pcl::PointXYZ pnt3 = cloud->points[p3];

				// get vectors vct1, vct2
				pcl::PointXYZ vct1;
				pcl::PointXYZ vct2;
				vct1.x = pnt2.x - pnt1.x;
				vct1.y = pnt2.y - pnt1.y;
				vct1.z = pnt2.z - pnt1.z;

				vct2.x = pnt3.x - pnt1.x;
				vct2.y = pnt3.y - pnt1.y;
				vct2.z = pnt3.z - pnt1.z;

				// calculate cross product v1 * v2 = normal vector nrm_vct
				pcl::PointXYZ nrm_vct;
				nrm_vct.x = (vct1.y * vct2.z) - (vct1.z * vct2.y);
				nrm_vct.y = (vct1.z * vct2.x) - (vct1.x * vct2.z);
				nrm_vct.z = (vct1.x * vct2.y) - (vct1.y * vct2.x);

				// calculate the plane equation coeffiecients
				qtn.pln_cff_a = nrm_vct.x;
				qtn.pln_cff_b = nrm_vct.y;
				qtn.pln_cff_c = nrm_vct.z;
				qtn.pln_cff_d = - ( (nrm_vct.x * (vct1.z - vct2.y)) +
				 					(nrm_vct.y * (vct1.x - vct2.z)) +
				 					(nrm_vct.z * (vct1.y - vct2.x)) );

				// calc the distance from plane equation to every other point
				// and store it in unordered_set if the distanceTol will fit
				for( int idx_dist = 0; idx_dist < cloud->points.size(); idx_dist++ )
				{
					pcl::PointXYZ pnt_dist = cloud->points[idx_dist];
					dist_cnt = 	qtn.pln_cff_a * pnt_dist.x + 
								qtn.pln_cff_b * pnt_dist.y + 
								qtn.pln_cff_c * pnt_dist.z +
								qtn.pln_cff_d;
					if( dist_cnt < 0 )
						dist_cnt *= -1;
					dist_sqrt = sqrt( qtn.pln_cff_a * qtn.pln_cff_a + 
									  qtn.pln_cff_b * qtn.pln_cff_b +
									  qtn.pln_cff_c * qtn.pln_cff_c );
					dist = dist_cnt / dist_sqrt;

					// evaluate distance tolerance
					if( dist < distanceTol )
					{
						inliersResult[inlier_idx].insert( idx_dist );
					}
				}
				if( inlier_idx == 1 )
				{
					if( inliersResult[1].size() > inliersResult[0].size() )
					{
						inliersResult[0].clear();
						inlier_idx = 0;
						inlier_out = 1;
					}
					else
					{
						inliersResult[1].clear();
						inlier_idx = 1;
						inlier_out = 0;
					}
				}
				else
				{
					if( inliersResult[0].size() > inliersResult[1].size() )
					{
						inliersResult[1].clear();
						inlier_idx = 1;
						inlier_out = 0;
					}
					else
					{
						inliersResult[0].clear();
						inlier_idx = 0;
						inlier_out = 1;
					}
				}
	 	//	}
	// 	}
	}

	return inliersResult[inlier_out];
}


std::unordered_set<int> Ransac(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, int maxIterations, float distanceTol )
{
	std::unordered_set<int> inliersResult[2];
	srand(time(NULL));
	double dist, dist_cnt, dist_sqrt;
	
	// TODO: Fill in this function
	struct
	{
		float ln_cff_a;
		float ln_cff_b;
		float ln_cff_c;
	} qtn;

	int inlier_result_idx = 0;
	int inlier_result_out = 0;


#ifdef NO
	for( idx_one = 0; idx_one < cloud->points.size(); idx_one++ )
	{
		// Combination without repetition = 190
		idx_two = idx_one + 1;
		for( ; idx_two < cloud->points.size(); idx_two++ )
		{
			// get the line equation
			pcl::PointXYZ pnt_one = cloud->points[idx_one];
			pcl::PointXYZ pnt_two = cloud->points[idx_two];
			qtn.ln_cff_a = pnt_one.y - pnt_two.y;
			qtn.ln_cff_b = pnt_two.x - pnt_one.x;
			qtn.ln_cff_c = pnt_one.x * pnt_two.y - pnt_two.x * pnt_one.y;
		
			// calc the distance from line equation to every other point
			// and store it in unordered_set if the distanceTol will fit
			for( int idx_dist = 0; idx_dist < cloud->points.size(); idx_dist++ )
			{
				//if( idx_dist != idx_one && idx_dist != idx_two )
				{
					pcl::PointXYZ pnt_dist = cloud->points[idx_dist];
					dist_cnt = qtn.ln_cff_a * pnt_dist.x + qtn.ln_cff_b * pnt_dist.y + qtn.ln_cff_c;
					if( dist_cnt < 0 )
						dist_cnt *= -1;
					dist_sqrt = sqrt( qtn.ln_cff_a * qtn.ln_cff_a + qtn.ln_cff_b * qtn.ln_cff_b );
					dist = dist_cnt / dist_sqrt;

					//std::cout << qtn.ln_cff_a << "\t" << qtn.ln_cff_b << "\t" << qtn.ln_cff_c << "\t" << 
					//			 "\t" << dist_cnt << "\t" << "\t" << dist_sqrt << "\t" << "\t" << dist << std::endl;
			
					// prove of distance tolerance
					if( dist < distanceTol )
					{
						inliersResult[inlier_result_idx].insert( idx_dist );
					}
				}
			}
			if( inlier_result_idx == 1 )
			{
				if( inliersResult[1].size() > inliersResult[0].size() )
				{
					inliersResult[0].clear();
					inlier_result_idx = 0;
					inlier_result_out = 1;
				}
				else
				{
					inliersResult[1].clear();
					inlier_result_idx = 1;
					inlier_result_out = 0;
				}
			}
			else
			{
				if( inliersResult[0].size() > inliersResult[1].size() )
				{
					inliersResult[1].clear();
					inlier_result_idx = 1;
					inlier_result_out = 0;
				}
				else
				{
					inliersResult[0].clear();
					inlier_result_idx = 0;
					inlier_result_out = 1;
				}
			}
		}
	}
#endif

	int iter = maxIterations;
	int pnt_rand_one, pnt_rand_two;
	while( --iter >= 0  )
	{
		pnt_rand_one = rand() % 20;
		pnt_rand_two = rand() % 20;
		if( pnt_rand_one == pnt_rand_two ) 
			continue;

		// get the line equation
		pcl::PointXYZ pnt_one = cloud->points[pnt_rand_one];
		pcl::PointXYZ pnt_two = cloud->points[pnt_rand_two];
		qtn.ln_cff_a = pnt_one.y - pnt_two.y;
		qtn.ln_cff_b = pnt_two.x - pnt_one.x;
		qtn.ln_cff_c = pnt_one.x * pnt_two.y - pnt_two.x * pnt_one.y;
		
		// calc the distance from line equation to every other point
		// and store it in unordered_set if the distanceTol will fit
		for( int idx_dist = 0; idx_dist < cloud->points.size(); idx_dist++ )
		{
			pcl::PointXYZ pnt_dist = cloud->points[idx_dist];
			dist_cnt = qtn.ln_cff_a * pnt_dist.x + qtn.ln_cff_b * pnt_dist.y + qtn.ln_cff_c;
			if( dist_cnt < 0 )
				dist_cnt *= -1;
			dist_sqrt = sqrt( qtn.ln_cff_a * qtn.ln_cff_a + qtn.ln_cff_b * qtn.ln_cff_b );
			dist = dist_cnt / dist_sqrt;

			// prove of distance tolerance
			if( dist < distanceTol )
			{
				inliersResult[inlier_result_idx].insert( idx_dist );
			}
		}

		if( inlier_result_idx == 1 )
		{
			if( inliersResult[1].size() > inliersResult[0].size() )
			{
				inliersResult[0].clear();
				inlier_result_idx = 0;
				inlier_result_out = 1;
			}
			else
			{
				inliersResult[1].clear();
				inlier_result_idx = 1;
				inlier_result_out = 0;
			}
		}
		else
		{
			if( inliersResult[0].size() > inliersResult[1].size() )
			{
				inliersResult[1].clear();
				inlier_result_idx = 1;
				inlier_result_out = 0;
			}
			else
			{
				inliersResult[0].clear();
				inlier_result_idx = 0;
				inlier_result_out = 1;
			}
		}
	}


	// Randomly sample subset and fit line

	// Measure distance between every point and fitted line
	// If distance is smaller than threshold count it as inlier

	// Return indicies of inliers from fitted line with most inliers
	
	return inliersResult[inlier_result_out];
}


int main ()
{

	// Create viewer
	pcl::visualization::PCLVisualizer::Ptr viewer = initScene();

	// Create data
	//pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = CreateData();
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = CreateData3D();
	

	// TODO: Change the max iteration and distance tolerance arguments for Ransac function
	//std::unordered_set<int> inliers = Ransac(cloud, 50, 0.5);
	std::unordered_set<int> inliers = RansacPlane(cloud, 100, 0.5);

	pcl::PointCloud<pcl::PointXYZ>::Ptr cloudInliers(new pcl::PointCloud<pcl::PointXYZ>());
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloudOutliers(new pcl::PointCloud<pcl::PointXYZ>());

	for(int index = 0; index < cloud->points.size(); index++)
	{
		pcl::PointXYZ point = cloud->points[index];
		//std::cerr << "x" << index << ": " << point.x << "\t" << "\t" << "y" << index << ": " << point.y << std::endl;

		if(inliers.count(index))
			cloudInliers->points.push_back(point);
		else
			cloudOutliers->points.push_back(point);
	}


	// Render 2D point cloud with inliers and outliers
	if(inliers.size())
	{
		renderPointCloud(viewer,cloudInliers,"inliers",Color(0,1,0));
		renderPointCloud(viewer,cloudOutliers,"outliers",Color(1,0,0));
	}
  	else
  	{
		renderPointCloud(viewer,cloud,"data");
  	}
	
  	while (!viewer->wasStopped ())
  	{
  	  viewer->spinOnce ();
  	}
}
