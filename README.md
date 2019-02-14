# BoundaryLayerLES

## Destabilizing buoyancy flux

Setups:

1. With a constant, destabilizing buoyancy flux:
- Initial linear stratification (`/destabilizing_buoyancy_flux/linear_stratification`)
- Initial mixed layer + linear stratification below (`/destabilizing_buoyancy_flux/mixed_layer`)
  * Test both AMD and Constant Smagorinsky

2. With a daily cycle of destabilizing buoyancy flux:
- With no internal heating during the day (convection shutoff)
- With internal heating during the day (shallow restratification)
