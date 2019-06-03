#ifndef VISALINGHAM_WHYATT_HPP
#define VISALINGHAM_WHYATT_HPP


#include <vector>

class VWSimplify {
    
    public:
        VWSimplify(const std::vector<std::pair<double,double>>& coords_in);
        
        const std::vector<size_t>& ordered_weights() const { return ordered_weights_; }
        const std::vector<std::tuple<double,double,double>>& coords() const { return coords_; }
        
        std::vector<std::pair<double,double>> by_threshold(double thres);     
        
    
        
    
    private:
        std::vector<std::tuple<double,double,double>> coords_;
        std::vector<size_t> ordered_weights_;
        
};
    


#endif
