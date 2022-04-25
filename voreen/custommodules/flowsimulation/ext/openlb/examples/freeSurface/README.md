# Free Surface  

To add a free surface behaviour to the program two general steps must be taken.  
The simple part is setting up the post processor.  
This example shows how to setup the first step in a 2D example.  

```
	FreeSurface2DSetup<T,DESCRIPTOR> free_surface_setup {
		superLattice, // Expecting a SuperLattice<T,DESCRIPTOR> here
		FreeSurface2D::Variables<T> {
			drop_isolated_cells, // Whether to convert interface cells to gas with no fluid and no interface cells neighbours
			transitionThreshold, // The anti jitter value when converting interface cells to gas or liquid
			lonelyThreshold, // When to convert interface cells which are lonely (No fluid neighbours or no gas neighbours)
			// Setting it to 0 always converts lonely cells. 1 would already be a value where a regular simulation wouldn't convert lonely cells
			has_surface_tension, // If surface tension should be active
			lattice_surface_tension_coefficient, // Surface tension coefficient in lattice units  
			force_conversion_factor, // Conversion value for the massless force (Don't use the converter function)  
			physical_delta // distance between two neighbouring cells in SI units on one axis  
		}
	};

	free_surface_setup.addPostProcessor();
```

The second more complicated setups occurs when creating the initial mass, epsilon and cell types.  
Due to the initialisation structure of openlb we need the physical distance in one axis direction, since neighbouring cells are not accessed easily.  
Generally the walls and borders are setup as usually. The fluid simulation area by default should be set to the gas type. Correspondingly the mass
and epsilon are set to zero.  
The seperation between the fluid and the gas phase is accomplished by a one layer thick interface layer.  
So after checking the rule whether to set the cell as fluid, every cell needs to check if this rule applies to their neighbour as well after not matching
the fluid rule. Then this cell is set as an interface cell. Otherwise the cells stays set as a gas cell.  
For the fluid cells the epsilon should be set to the value 1.    
The simplest way to achieve this is to build an analytical functor. The 2D falling drop example is set as follows.  

```
template <typename T, typename DESCRIPTOR>
class FreeSurfaceFallingDrop2D final : public AnalyticalF2D<T,T> {
private:
  T lattice_size;
  std::array<T, 3> cell_values;
public:
  FreeSurfaceFallingDrop2D(T lattice_size_, const std::array<T,3>& cell_vals):AnalyticalF2D<T,T>{1}, lattice_size{lattice_size_}, cell_values{cell_vals}{}

  bool operator()(T output[], const T x[]) override {
    output[0] = cell_values[0];
    T radius = 0.00155;

    if(x[1] <= radius){
      output[0] = cell_values[2];
    }else if(x[1] <= radius + lattice_size * 1.1){
      output[0] = cell_values[1];
    }

    std::array<T,DESCRIPTOR::d> point = {0.015, 2 * radius + lattice_size * 4};
    std::array<T,DESCRIPTOR::d> diff = {x[0] - point[0], x[1] - point[1]};

    if((diff[0]*diff[0] + diff[1] * diff[1]) <= radius*radius){
      output[0] = cell_values[2];
    }else{
      for(int i = -1; i <= 1; ++i){
        for(int j = -1; j <= 1; ++j){
          std::array<T,DESCRIPTOR::d> shifted_diff = {diff[0]+i*lattice_size*1.1, diff[1]+j*lattice_size*1.1};
          if((shifted_diff[0]*shifted_diff[0] + shifted_diff[1] * shifted_diff[1]) <= radius*radius){
            output[0] = cell_values[1];
            return true;
          }
        }
      }
    }

    return true;
  }
};
```

This is later used in the lattice preparation step as 

```
void prepareFallingDrop(UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice<T, DESCRIPTOR>& sLattice,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     SuperGeometry<T,2>& superGeometry, T lattice_size, const FreeSurfaceAppHelper& helper)
{
  AnalyticalConst2D<T,T> zero( 0. );
  AnalyticalConst2D<T,T> one( 1. );
  AnalyticalConst2D<T,T> two( 2. );
  AnalyticalConst2D<T,T> four( 4. );
  FreeSurfaceFallingDrop2D<T,DESCRIPTOR> cells_analytical{ lattice_size, {0., 1., 2.}};
  FreeSurfaceFallingDrop2D<T,DESCRIPTOR> mass_analytical{ lattice_size, {0., 0.5, 1.}};

  AnalyticalConst2D<T,T> force_zero{0., 0.};

  for (int i: {0,1,2}) {
    sLattice.defineField<FreeSurface::MASS>(superGeometry, i, zero);
    sLattice.defineField<FreeSurface::EPSILON>(superGeometry, i, zero);
    sLattice.defineField<FreeSurface::CELL_TYPE>(superGeometry, i, zero);
    sLattice.defineField<FreeSurface::CELL_FLAGS>(superGeometry, i, zero);
    sLattice.defineField<descriptors::FORCE>(superGeometry, i, force_zero);
  }

  sLattice.defineField<FreeSurface::CELL_TYPE>(superGeometry, 1, cells_analytical);

  sLattice.defineField<FreeSurface::MASS>(superGeometry, 1, mass_analytical);
  sLattice.defineField<FreeSurface::EPSILON>(superGeometry, 1, mass_analytical);

  for (int i: {0,2}) {
    sLattice.defineField<FreeSurface::EPSILON>(superGeometry, i, one);
    sLattice.defineField<FreeSurface::CELL_TYPE>(superGeometry, i, four);
  }

  T force_factor = 1./ converter.getConversionFactorForce() * converter.getConversionFactorMass();
  AnalyticalConst2D<T,T> force_a{helper.gravity_force[0] * force_factor, helper.gravity_force[1] * force_factor};
  sLattice.defineField<descriptors::FORCE>(superGeometry.getMaterialIndicator({1}), force_a);

}
```

Here the gravity force is set as well. These steps are sufficient to setup a basic free surface simulation though additional steps can be taken.  
The Smagorinsky BGK dynamics classes `SmagorinskyBGKdynamics<T,DESCRIPTOR>` or `SmagorinskyForcedBGKdynamics<T,DESCRIPTOR>` are normally used for
free surface setups.  
Additionally the epsilon values of walls are set to 1 in the falling drop example. This influences if the wall is wetting or not. With 1 being a wetting wall.  


