#ifndef _HS_SFM_SFM_FILE_IO_TRACKS_SAVER_HPP_
#define _HS_SFM_SFM_FILE_IO_TRACKS_SAVER_HPP_

#include <string>
#include <fstream>

#include "hs_sfm/sfm_utility/match_type.hpp"

namespace hs
{
namespace sfm
{
namespace fileio
{

struct TracksSaver
{
  typedef hs::sfm::TrackContainer TrackContainer;
  typedef hs::sfm::ViewInfo ViewInfo;
  typedef hs::sfm::ViewInfoIndexer ViewInfoIndexer;
  typedef int Err;
  Err operator() (const std::string& tracks_path,
                  const TrackContainer& tracks,
                  const ViewInfoIndexer& view_info_indexer) const
  {
    std::ofstream tracks_file(tracks_path.c_str(), std::ios::out);
    if (!tracks_file.is_open())
    {
      return -1;
    }

    size_t number_of_tracks = tracks.size();
    tracks_file<<number_of_tracks<<"\n";
    for (size_t track_id = 0; track_id < number_of_tracks; track_id++)
    {
      size_t number_of_views = tracks[track_id].size();
      tracks_file<<number_of_views;
      for (size_t i = 0; i < number_of_views; i++)
      {
        size_t image_id = tracks[track_id][i].first;
        size_t key_id = tracks[track_id][i].second;
        const ViewInfo* view_info =
          view_info_indexer.GetViewInfoByTrackImage(track_id, image_id);
        int flag = view_info->is_blunder ? 1 : 0;
        tracks_file<<" "<<image_id<<" "<<key_id<<" "<<flag;
      }
      tracks_file<<"\n";
    }

    return 0;
  }
};

}
}
}

#endif
