/*
 * Copyright (C) 2019 - 2020 by the emerging fields initiative 'Novel Biopolymer
 * Hydrogels for Understanding Complex Soft Tissue Biomechanics' of the FAU
 *
 * This file is part of the EFI library.
 *
 * The EFI library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * Author: Stefan Kaessmair
 */

// efi headers
#include <efi/worker/boundary_worker.h>
#include <efi/factory/registry.h>


using namespace dealii;

namespace efi
{

template <int dim>
void
BoundaryWorker<dim>::
do_fill (ScratchData<dim> & scratch_data,
         CopyData         & copy_data) const
{
    using namespace dealii;
    using namespace dealii::Physics::Elasticity;

    // std::cout << "Entered boundary loop" << std::endl;

    auto global_vector_name = Extractor<dim>::global_vector_name();

    // Get the number of quadrature points and
    // the number of dofs per cell.
    auto q_points    = ScratchDataTools::get_quadrature_points (scratch_data); // 4 q_points
    

    auto n_q_points    = ScratchDataTools::n_quadrature_points (scratch_data);
    // auto dofs_per_cell = ScratchDataTools::dofs_per_cell (scratch_data);

    // Get the current fe values object, which
    // is needed to access the shape functions.
    auto &fe  = static_cast<const dealii::FEFaceValuesBase<dim>&> (
             ScratchDataTools::get_current_fe_values (scratch_data));
    auto &JxW = ScratchDataTools::get_JxW_values        (scratch_data);
    auto &N   = ScratchDataTools::get_normal_vectors (scratch_data);
    auto &load   = ScratchDataTools::get_load (scratch_data);
    auto dofs_per_cell = fe.dofs_per_cell;

    
    // auto &Grad_u = ScratchDataTools::get_gradients (
    //                        scratch_data, global_vector_name,
    //                        Extractor<dim>::displacement(),ad_type(0));
    // auto &F   = ScratchDataTools::get_deformation_grads  (scratch_data, global_vector_name, ad_type());
    // auto &J   = ScratchDataTools::get_jacobians (scratch_data,global_vector_name,ad_type());

    // Create references to the copy data objects
    // for easy accessibility.
    auto &local_rhs    = copy_data.vectors.back();
    // auto &local_matrix = copy_data.matrices.back();

    local_rhs = 0;
    // local_matrix = 0;

    // loop over all quadrature points
    for (unsigned int q = 0; q < n_q_points; ++q)
    {
        dealii::Tensor<1, dim> traction = -1*load*(N[q]);
        
        // loop over all degrees of freedom (rows)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
            const unsigned int component_i = fe.get_fe().system_to_component_index(i).first;
            double Ni = fe.shape_value(i,q);
            local_rhs(i) += Ni*traction[component_i]*JxW[q];  
        }// i loop
    }// q loop
}


// Instantiation
template class BoundaryWorker<2>;
template class BoundaryWorker<3>;

// Registration
EFI_REGISTER_OBJECT (EFI_TEMPLATE_CLASS (BoundaryWorker,2));
EFI_REGISTER_OBJECT (EFI_TEMPLATE_CLASS (BoundaryWorker,3));

}//namespace efi


