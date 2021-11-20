for (std::size_t i = 0; i < clusters.size(); ++ i)
    {
      if(clusters[i].is_free)
        {
          clusters[i].is_free = false;
          FT max_area = clusters[i].area;
          std::size_t index_max_area = i;

          //initialization containers
          std::vector < std::size_t > index_container;
          index_container.push_back(i);
          std::vector < std::size_t > index_container_former_ring;
          index_container_former_ring.push_back(i);
          std::list < std::size_t > index_container_current_ring;

          //propagation
          bool propagation=true;
          do
            {
              propagation=false;

              //neighbors
              for (std::size_t k=0;k<index_container_former_ring.size();k++)
                {

                  std::size_t cluster_index=index_container_former_ring[k];

                  for (std::size_t j = 0; j < clusters[cluster_index].orthogonal_clusters.size(); ++ j)
                    {
                      std::size_t cluster_index_2 = clusters[cluster_index].orthogonal_clusters[j];
                      if(clusters[cluster_index_2].is_free)
                        {
                          propagation = true;
                          index_container_current_ring.push_back(cluster_index_2);
                          clusters[cluster_index_2].is_free = false;

                          if(max_area < clusters[cluster_index_2].area)
                            {
                              max_area = clusters[cluster_index_2].area;
                              index_max_area = cluster_index_2;
                            }
                        } 
                    }
                }

              //update containers
              index_container_former_ring.clear();
              for(std::list < std::size_t>::iterator it = index_container_current_ring.begin();
                  it != index_container_current_ring.end(); ++it)
                {
                  index_container_former_ring.push_back(*it);
                  index_container.push_back(*it);
                }
              index_container_current_ring.clear();

            }
          while(propagation);