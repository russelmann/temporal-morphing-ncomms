#ifndef EIGEN_DENSE_BASE_ADDONS_H
#define EIGEN_DENSE_BASE_ADDONS_H

// This is included only with cereal; change if other addons needed

template<class Archive>
void save(Archive & ar) const
{
	derived().eval();
	ar( cereal::make_nvp("rows", derived().rows()) );
	ar( cereal::make_nvp("cols", derived().cols()) );
	ar( cereal::make_nvp("data", std::vector<Scalar>(derived().data(), derived().data() + derived().size())) );
}

template<class Archive>
void load(Archive & ar)
{
	Index rows;
	Index cols;
	ar( cereal::make_nvp("rows", rows) );
	ar( cereal::make_nvp("cols", cols) );
	if (rows != derived().rows() || cols != derived().cols())
		derived().resize(rows, cols);
	
	if (rows == 0 || cols == 0) return;

	std::vector<Scalar> vec;
	ar( cereal::make_nvp("data", vec) );
	memcpy(derived().data(), vec.data(), vec.size() * sizeof(Scalar));
}

#endif