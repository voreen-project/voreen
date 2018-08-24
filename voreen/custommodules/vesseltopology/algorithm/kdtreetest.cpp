#define TGT_ASSERT_H
#define tgtAssert(cond, text) assert((cond))
#include "tgt/vector.h"
#include "kdtree.h"

#include <random>
#include <iostream>


struct ExampleElement {
    typedef float CoordType;
    const tgt::Vector3<CoordType>& getPos() const {
        return pos_;
    }

    ExampleElement(tgt::vec3 pos, int data)
        : pos_(pos)
        , data_(data)
    {
    }
    tgt::vec3 pos_;
    int data_;
};

int main() {

    std::default_random_engine generator;
    std::uniform_real_distribution<float> pos_distribution(0,1);

    voreen::kdtree::ElementArrayBuilder<ExampleElement> builder("elms.tmp");
    for(int i=0; i<10000; ++i) {
        builder.push(ExampleElement(tgt::vec3(
                        pos_distribution(generator),
                        pos_distribution(generator),
                        pos_distribution(generator)),
                    i));
    }

    voreen::kdtree::Tree<ExampleElement> tree("tree.bin", std::move(builder));
    std::cout << tree.root().elm_.data_ << std::endl;

}
