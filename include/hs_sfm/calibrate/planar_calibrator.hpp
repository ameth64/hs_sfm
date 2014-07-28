#ifndef _HS_SFM_CALIBRATE_PLANAR_CALIBRATOR_HPP_
#define _HS_SFM_CALIBRATE_PLANAR_CALIBRATOR_HPP_
#include"opencv2\calib3d.hpp"
#include"opencv2\core.hpp"
#include"opencv2\features2d.hpp"
#include"opencv2\highgui.hpp"
#include"opencv2\imgproc.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"

namespace hs
{
namespace sfm
{
namespace calibrate
{

template <typename _Scalar>
class PlanarCalibrator
{
public:
  typedef _Scalar Scalar;
  typedef hs::sfm::CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef hs::sfm::CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef EIGEN_MATRIX(Scalar, 2, 2) KeyCovariance;
  typedef EIGEN_MATRIX(Scalar, 3, 3) PointCovariance;
  typedef EIGEN_VECTOR(Scalar, 2) Key;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef std::pair<Key, Point> Correspondence;
  typedef EIGEN_STD_VECTOR(Correspondence) PatternView;
  typedef EIGEN_STD_VECTOR(PatternView) PatternViewContainer;
  typedef size_t ImageWidth;
  typedef size_t ImageHigth;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef hs::math::geometry::Rotation3D<Scalar> Rotation;

public:
	int operator() (const PatternViewContainer& pattern_views,
				  const ImageWidth &image_width,
	              const ImageHigth &image_higth,
                  const KeyCovariance& key_covariance,
				  const PointCovariance& point_covariance,
                  IntrinsicParams& intrinsic_params,
                  ExtrinsicParamsContainer& extrinsic_params_set) const
  {
		std::vector<cv::Point3f> realPoint;
		std::vector<cv::Point2f> cornerPoint;
		std::vector<std::vector<cv::Point3f> > objectPoint;
		std::vector<std::vector<cv::Point2f> > imagePoint;
		for(int image_count=0;image_count!=pattern_views.size();image_count++)
		{
			for(int point_count=0;point_count!=pattern_views[image_count].size();point_count++)
			{
				realPoint.push_back(cv::Point3f((pattern_views[image_count][point_count]).second(0),pattern_views[image_count][point_count].second(1),pattern_views[image_count][point_count].second(2)));
				cornerPoint.push_back(cv::Point2f(pattern_views[image_count][point_count].first(0),pattern_views[image_count][point_count].first(1)));
			}
			objectPoint.push_back(realPoint);
			imagePoint.push_back(cornerPoint);
		}
		std::vector<cv::Mat> tvecs,rvecs;
		cv::Mat cameraMatrix=cv::Mat::eye(3,3,CV_64F);
		cv::Mat distCoeffs= cv::Mat::zeros(1, 5, CV_64F);
		cv::calibrateCamera(objectPoint,imagePoint,cv::Size(image_width,image_higth),cameraMatrix,distCoeffs,rvecs,tvecs);

		intrinsic_params = IntrinsicParams(Scalar(cameraMatrix.at<float>(0, 0)), 
			Scalar(0), Scalar(cameraMatrix.at<Scalar>(0, 2)),
			(Scalar)cameraMatrix.at<Scalar>(1, 2), (Scalar)1,
			(Scalar)distCoeffs.at<Scalar>(0, 0),
			(Scalar)distCoeffs.at<Scalar>(0, 1),
			(Scalar)distCoeffs.at<Scalar>(0, 4),
			(Scalar)distCoeffs.at<Scalar>(0, 2),
			(Scalar)distCoeffs.at<Scalar>(0, 3)
	);

		for (int i = 0; i != tvecs.size(); i++)
		{

			extrinsic_params_set.push_back(ExtrinsicParams
														(Rotation
															(Vector3
																((Scalar)rvecs[i].at<Scalar>(0, 0),
																(Scalar)rvecs[i].at<Scalar>(1, 0),
																(Scalar)rvecs[i].at<Scalar>(2, 0)
															    )
															),
															 Vector3
																(tvecs[i].at<Scalar>(0, 0),
																tvecs[i].at<Scalar>(1, 0),
																tvecs[i].at<Scalar>(2, 0)
																 )
														)
				);
		}
    return 0;
  }
};

}
}
}

#endif
