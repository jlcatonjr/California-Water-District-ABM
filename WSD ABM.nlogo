
; WATER DISTRICT ABM BY MEGA P

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; EXTENSIONS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

extensions [matrix GIS bitmap]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; AGENT BREEDS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; breed [farmers farmer]
breed [drillers driller]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; GLOBAL VARIABLES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

globals
[ 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;GROUNDWATER MODEL;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
  max-error   
  max-head
  min-head
  max-K
  min-K
  iterations
  aquifer-width
  A
  A-row-list
  A-row-reduced-list
  A-row-reduced-matrix
  A-column-list
  A-column-reduced-list
  A-reduced
  inverse-A
  C
  C-row-list
  C-row-reduced-list
  C-row-reduced-matrix
  C-reduced
  solution-vector
  number-of-unknowns
  K-patch-read

  ;;;;;;;;lists of patches;;;;;;;
  no-flow-cells-list
  fixed-head-cells-list
  remove-cells-list
  remove-cells-list-corrected
  active-cells-remap
  active-cells-remap-indexes
  
  ;;;;;;;multiplier functions;;;;;;;
  sine-recharge
  sine-well
  
  ;;;;;;river budget term;;;;;;;;
  riv-budget
  
  ;;;;;;well pumping budget term;;;;;;;;
  well-budget

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;SOCIAL MODEL;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;
  
  year
  month

  base-water-price
  variable-water-price
  water-price 

  canal-water-used
  well-water-used

  wd-artichoke 
  wd-broccoli 
  wd-cabbage 
  partichoke 
  pbroccoli 
  pcabbage 
  yartichoke 
  ybroccoli 
  ycabbage 
  ocartichoke 
  ocbroccoli 
  occabbage 
  crop-list
  profit-list 
  yield-list 
  water-price-list 
  wd-list 
  oc-list
  waiting-period-list
  count-bankrupt
  
  cropnames-list
  acreage-list
  acreage-list-pct
  wdAF-list
  wdM3-list
  netprice-list
  
  current-canal-volume
  canal-water-price
  
  max-depth-list           ;; list containing maximum water table depths historically observed in the aquifer
  max-depth-value          ;; the maximum historical water table depth recorded (this is used by the driller to determine depths on new wells)

;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;Crop Samples;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;
    
  ;;;;;;;CROP AREAS;;;;;;;;;;
  
  %-of-Grains 
  %-of-Alfalfa 
  %-of-Cotton 
  %-of-Onion 
  %-of-Truck 
  %-of-Deciduous 
  %-of-Grapes 
  %-of-Potatoes 
  %-of-Melon 
  %-of-Carrots 
  %-of-Citrus 
  %-of-Fallow
 
  output-almonds
  output-cherries
  output-citrus
  output-grapes
  output-pistachios
  output-tomatoes
  output-cotton
  output-alfalfa
  output-silageandforage
  output-onions
  output-bellpeppers
  output-potatoes
  output-fallow
  
  output-almonds_name
  output-cherries_name
  output-citrus_name
  output-grapes_name
  output-pistachios_name
  output-tomatoes_name
  output-cotton_name
  output-alfalfa_name
  output-silageandforagee
  output-onions_name
  output-bellpeppers_name
  output-potatoes_name
  output-fallow_name
  
  filename
  cost-differentials_name
  price-sensitivity_name 
  mutation-rate_name
  starting-water-price_name
  electricity-price_name 
  %-optimizers_name
  target-water-level_name
  canal-water-price_name 
  drought-index_name 
  prices_name
  inflow-to-basin_name 

;  establishment-revenue-list
;  establishment-costs-list
;  establishment-water-demand-list

  ready-to-histogram?
  m-list
  n-list
  
  total-wells-drilled
  
  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;DATASETS;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;
  
  csv
  fileList
  
  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;GINI;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;
  
  gini-index-reserve
  lorenz-points  
  model-ready?
  plots-ready?
  
  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;; GIS;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;
  
  kern-shapefile  
  kern-raster
  gis-envelope
  
]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; AGENT/PATCH VARIABLES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

patches-own
[ 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;GROUNDWATER MODEL;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ;;;;;;;;source and sink terms;;;;;;;
 
  Qwell                                    ;; well pumping rate                                        [L3/T] 
  Qwell-temp                               ;; temporary dummy variable to initialize heads
  Qcanal
  Qcanal-temp
  
  
  Qinjection                               ;; well injection rate                                      [L3/T] 
  Qinjection-temp                          ;; temporary dummy variable to initialize heads
 
  Recharge                                 ;; recharge                                                 [L3/L2*T] ~ [L/T] 
  Recharge-temp                            ;; temporary dummy variable to initialize heads
  
  ET                                       ;; evapotranspiration, function of depth to water table     [L3/L2*T] ~ [L/T] 
  ET-max-rate-temp                         ;; temporary dummy variable to initialize heads
  ET-max-rate-patch
  ET-extinction-depth-patch
  ET-land-surface-elevation-patch
  
  DRAIN                                    ;; drain discharge, depends on h-d             [L3/T] Note: divide by cell area in sourceterm equation
  DRAIN-elevation-patch 
  DRAIN-conductance-patch 
  
  RIV                                      ;; leakage from river                          [L3/L2*T] ~ [L/T] 
  RIV-elevation-patch
  RIV-bottom-patch
  RIV-conductance-patch
  
  ;;;;;;patch hydraulic parameters;;;;;;;
  S-patch                                  ;; storativity of each patch
  S-patch-temp                             ;; temporary dummy variable to solve for steady state
  K-patch                                  ;; conductivity of each patch
  T-patch                                  ;; transmissivity of each patch

  ;;;;;;;;;type of patch;;;;;;;;
  interior-node?                           ;; is this an active patch?
  adjacent-to-boundary?                    ;; is this patch adjacent to a boundary?
  interior-but-corner?                     ;; out of the interior nodes, is it in the corner?
  
  ;;;;;;patch bc's;;;;;
  no-flow?                                 ;; is this a no-flow cell?
  fixed-head?                              ;; is this a fixed head cell?
  fixed-flux?                              ;; is this a fixed flux cell?
  
  ;;;;;;patch neighbours;;;;
  N-neighbour
  S-neighbour
  E-neighbour
  W-neighbour
  
  ;;;;;;patch features;;;;;
  well?
  injection?
  recharge?
  ET?
  DRAIN?
  RIV?
  
  ;;;;;;;head values;;;;;;
  H                                       ;; heads at the new timestep
  Hinew                                   ;; heads at the new iteration
  Hinitial                                ;; a reference H to calculate drawdowns
  Hground                                 ;; level of ground
  H-ex-ante                               ;; level of ground water before planting
  ;;;;;;solution process parameters;;;;;
  TN
  TS
  TE
  TW
  TC
  Ncof
  Scof
  Ecof
  Wcof
  Ccof
  SourceTerm
  RHS
  patch-number-tag
  
  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;SOCIAL MODEL;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;
  
  water 
  money 
  wealth
  current-crop
  current-crop-index 
  
  gmlist
  gmlist-temp
  past-crop-index-list
  past-crop-list
  past-gm-list
  past-water-price-list
  weight-list
  wait-counter
  current-crop-temp
  position-of-this-crop
  mycrop
  bankrupt?
  water-level
  well-depth
  pumping-cost
  dry-well?
  new-well-needed?
  mortgage?
  payments-made
  
  canal-rights-rank
  canal-volume-claimed
  canal-water? 
  
  crop-assigned?
  farmer-well-here?
  
  gw-bank-here?
  farmer-here?
  optimizer?
  copycat?
  tree-grower?
  
  mutate?
]


drillers-own
[
  client-list                     ;;list of farmers waiting for a new well
  drilling-depth                  ;;how deep they drill — [m]
  drilling-duration               ;;how long they take to build a well — [months] 
  drilling-queue                  ;;how many farmers are waiting for a well
]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; GROUNDWATER MODEL VIEWS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to refresh-view
  if view = "patch numbering" [ask patches with [interior-node? = true] [set plabel precision patch-number-tag 1 set plabel-color black]]
  if view = "boundary conditions" [view-bc ask patches with [no-flow? = true][set pcolor black] ask patches with [fixed-head? = true][set pcolor blue]ask patches with [RIV? = true][set pcolor cyan] ask patches with [well? = true][set pcolor red] ask patches with [injection? = true][set pcolor blue]]
  if view = "wells" [ask patches with [interior-node? = true][set plabel ""] ask patches with [well? = true][set plabel-color black set plabel "W"] ask patches with [injection? = true][set plabel-color black set plabel "R"]]
  if view = "T (values)" [ask patches [set plabel "" set plabel precision T-patch 0 set plabel-color black]]
  if view = "K (values)" [ask patches [set plabel "" set plabel precision K-patch 0 set plabel-color black]]
  if view = "K (contours)" [setup-maxmin-K ask patches with [interior-node? = true and no-flow? = false and ET? = false and fixed-head? = false and DRAIN? = false and RIV? = false] [set pcolor scale-color brown K-patch min-K max-K set plabel precision K-patch 0 set plabel-color black] ask patches with [no-flow? = true][set pcolor black]]
  if view = "S (values)" [ask patches [set plabel "" set plabel precision S-patch 5 set plabel-color black]]

  if view = "crops" [ask patches [set plabel ""]]
  if view = "heads (values)" [ask patches [set plabel precision H 0 set plabel-color black]]
  if view = "heads (contours)" [setup-maxmin-heads ask patches with [interior-node? = true and no-flow? = false and ET? = false and fixed-head? = false and DRAIN? = false and RIV? = false] [setup-maxmin-heads set pcolor scale-color gray H min-head max-head set plabel precision H 0 set plabel-color black] ask patches with [no-flow? = true][set pcolor black] ask patches with [ET? = true][set plabel precision H 0 set plabel-color black]] 

  if view = "areal recharge" [ask patches [set plabel "" set plabel precision Recharge 4 set plabel-color black]]
  
  if view = "ET" [ask patches [set plabel "" ask patches with [ET? = true][set plabel precision Recharge 4 set plabel-color black set pcolor green]]]
  
  if view = "DRAIN conductance" [ask patches with [DRAIN? = true][set plabel DRAIN-conductance-patch]]
  if view = "DRAIN flow [m3/d] and head values" [ask patches with [DRAIN? = true][set plabel precision DRAIN 0] ask patches with [DRAIN? = false][set plabel precision H 0 set plabel-color black]]                    ;; cells that are not drains show the head value
  if view = "DRAIN flow [L/s] and head values" [ask patches with [DRAIN? = true][set plabel precision (DRAIN * 1000 / 86400) 0] ask patches with [DRAIN? = false][set plabel precision H 0 set plabel-color black]]    ;; cells that are not drains show the head value

  if view = "RIV flow [m3/d] and head values" [ask patches with [RIV? = true][set plabel precision RIV 0 if RIV < 0 [set pcolor red] if RIV >= 0 [set pcolor cyan]] ask patches with [RIV? = false][set plabel precision H 0 set plabel-color black]]                    ;; cells that are not drains show the head value
  if view = "RIV flow [L/s] and head values" [ask patches with [RIV? = true][set plabel precision (RIV * 1000 / 86400) 0 if RIV < 0 [set pcolor red] if RIV >= 0 [set pcolor cyan]] ask patches with [RIV? = false][set plabel precision H 0 set plabel-color black]]    ;; cells that are not drains show the head value
  if view = "RIV flow [L/s] and drawdowns" [ask patches with [RIV? = true][set plabel precision (RIV * 1000 / 86400) 0 if RIV < 0 [set pcolor red] if RIV >= 0 [set pcolor cyan]] ask patches with [RIV? = false][set plabel precision (Hinitial - H) 0 set plabel-color black]]    ;; cells that are not drains show the head value

  if view = "RIV stage and head values" [ask patches with [DRAIN? = true][set plabel precision RIV-elevation-patch 0] ask patches with [DRAIN? = false][set plabel precision H 0 set plabel-color black]]      ;; cells that are not drains show the head value
  if view = "RIV bottom and head values" [ask patches with [DRAIN? = true][set plabel precision RIV-bottom-patch 0] ask patches with [DRAIN? = false][set plabel precision H 0 set plabel-color black]]        ;; cells that are not drains show the head value
  if view = "RIV conductance" [ask patches with [DRAIN? = true][set plabel precision RIV-conductance-patch 0] ask patches with [DRAIN? = false][set plabel ""]]                                                ;; cells that are not drains are blank
end

to reset-view
  ask patches [set pcolor white set plabel "" set plabel-color black view-bc]
end


to view-bc
  ifelse left-bc = "no-flow" 
  [ask patches with [pxcor = 0][set pcolor black]]
  [ask patches with [pxcor = 0][set pcolor blue]]
  
  ifelse right-bc = "no-flow" 
  [ask patches with [pxcor = max-pxcor][set pcolor black]]
  [ask patches with [pxcor = max-pxcor][set pcolor blue]]
  
  ifelse top-bc = "no-flow" 
  [ask patches with [pycor = max-pycor][set pcolor black]]
  [ask patches with [pycor = max-pycor][set pcolor blue]]
  
  ifelse bottom-bc = "no-flow" 
  [ask patches with [pycor = min-pycor][set pcolor black]]
  [ask patches with [pycor = min-pycor][set pcolor blue]]
  
  ask patches with [fixed-head? = true][set pcolor blue]
  ask patches with [ET? = true][set pcolor brown]
  ask patches with [DRAIN? = true][set pcolor magenta]
  ask patches with [RIV? = true][set pcolor cyan]
end


to setup-maxmin-heads                                                    ;; finds max and min H values for colour-coding
  set max-head [H] of max-one-of patches with [interior-node? = true and no-flow? = false][H] 
  set min-head [H] of min-one-of patches with [interior-node? = true and no-flow? = false][H] 
end


to setup-maxmin-K                                                        ;; finds max and min K values for colour-coding
  set max-K [K-patch] of max-one-of patches [K-patch]
  set min-K [K-patch] of min-one-of patches [K-patch]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; SET INITIAL HEADS TO CALCULATE DRAWDOWNS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to set-initial-heads
  ask patches with [interior-node? = true][
    set Hground 500
    set Hinitial 400]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; SET BCs / POINT-AND-CLICK ON THE INTERFACE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to change-to-fixed-head
   if mouse-down?     ;; reports true or false to indicate whether mouse button is down
    [
      ask patch mouse-xcor mouse-ycor
      [ 
        set pcolor blue 
        set fixed-head? true
        set no-flow? false
        set H fixed-head-value-new
        set K-patch K 
        ifelse aquifer-type = "confined"
        [set T-patch (K-patch * aquifer-thickness)]       ;; sets transmissivity and storage for confined conditions
        [set T-patch (K-patch * H)]                      ;; sets transmissivity and storage for unconfined conditions   
        set plabel-color black
        set plabel H
      ]
    ]
end


to change-to-no-flow
   if mouse-down?     ;; reports true or false to indicate whether mouse button is down
    [
      ask patch mouse-xcor mouse-ycor
      [ 
        set pcolor black 
        set fixed-head? false
        set no-flow? true
      ]
    ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; SET HYDRAULIC PARAMETERS / POINT-AND-CLICK ON THE INTERFACE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to set-K
  if mouse-down?     ;; reports true or false to indicate whether mouse button is down
    [
      ask patch mouse-xcor mouse-ycor
      [ 
        set K-patch K-input 
        set plabel-color black
        set plabel K-patch
      ]
    ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; SET DRAIN / POINT-AND-CLICK ON THE INTERFACE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to set-DRAIN-patch
  if mouse-down?     ;; reports true or false to indicate whether mouse button is down
    [
      ask patch mouse-xcor mouse-ycor
      [ 
        set DRAIN? true
        set DRAIN-elevation-patch DRAIN-elevation
        set DRAIN-conductance-patch DRAIN-conductance
        set plabel-color black
        set pcolor magenta
        set plabel DRAIN-conductance-patch
      ]
    ]
end

to clear-DRAIN-patches-all
  ask patches with [DRAIN? = true]
  [
    set DRAIN? false 
    set DRAIN 0
    set DRAIN-elevation-patch 0
    set DRAIN-conductance-patch 0 
    set pcolor white
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; SET RIV / POINT-AND-CLICK ON THE INTERFACE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to set-RIV-patch
  if mouse-down?     ;; reports true or false to indicate whether mouse button is down
    [
      ask patch mouse-xcor mouse-ycor
      [ 
        set RIV? true
        set RIV-elevation-patch RIV-elevation
        set RIV-bottom-patch RIV-bottom
        set RIV-conductance-patch RIV-conductance
        set plabel-color black
        set pcolor cyan
        set plabel RIV-elevation-patch
      ]
    ]
end

to clear-RIV-patches-all
  ask patches with [RIV? = true]
  [
    set RIV? false 
    set RIV 0
    set RIV-elevation-patch 0
    set RIV-bottom-patch 0
    set RIV-conductance-patch 0
    set pcolor white
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; SET ET / POINT-AND-CLICK ON THE INTERFACE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to set-ET-patch
  if mouse-down?     ;; reports true or false to indicate whether mouse button is down
    [
      ask patch mouse-xcor mouse-ycor
      [ 
        set ET? true
        set ET-max-rate-patch ET-max-rate
        set ET-extinction-depth-patch ET-extinction-depth
        set ET-land-surface-elevation-patch ET-land-surface-elevation
        set plabel-color black
        set pcolor brown
        set plabel ET-max-rate-patch
      ]
    ]
end
 

to clear-ET-patches-all
  ask patches with [ET? = true]
  [
    set ET? false 
    set ET-max-rate 0
    set ET-extinction-depth 0
    set ET-land-surface-elevation 0
    set pcolor white
  ]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; SET FIXED FLUX — POINT-AND-CLICK ON THE INTERFACE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to place-fixed-flux
  if mouse-down?     ;; reports true or false to indicate whether mouse button is down
    [
      ask patch mouse-xcor mouse-ycor
        [ 
          set pcolor magenta 
          set fixed-flux? true
          set injection? true
          set Qinjection Qfixed-flux
          set plabel-color black
          set plabel Qinjection
        ]
    ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; SET WELLS — POINT-AND-CLICK ON THE INTERFACE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to place-pumping-well
  if mouse-down?     ;; reports true or false to indicate whether mouse button is down
    [
      ask patch mouse-xcor mouse-ycor
      [ 
        set pcolor red 
        set well? true
        set Qwell Qwell-input
        ;set screen-level screen-level-input
        set plabel-color black
        set plabel Qwell
      ]
    ]
end


to place-injection-well
  if mouse-down?     ;; reports true or false to indicate whether mouse button is down
    [
      ask patch mouse-xcor mouse-ycor
        [ 
          set pcolor blue 
          set injection? true
          set Qinjection Qinjection-input
          set plabel-color black
          set plabel Qinjection
        ]
    ]
end


;; the following routines allow clearing the view from stresses and resetting them using the interface buttons

to hide-stresses-patches
  ask patches with [well? = true or injection? = true]
  [
    if well? = true [set pcolor white set plabel-color red]
    if injection? = true [set pcolor white set plabel-color blue]
  ]
end


to hide-stresses-labels
  ask patches with [well? = true or injection? = true]
  [
    if well? = true [set plabel ""]
    if injection? = true [set plabel ""]
  ]
end


to clear-stresses
  ask patches with [well? = true or injection? = true]
  [
    if well? = true [set Qwell 0 hide-stresses-patches hide-stresses-labels set well? false]
    if injection? = true [set Qinjection 0 hide-stresses-patches hide-stresses-labels set injection? false]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; SET RECHARGE / POINT-AND-CLICK ON THE INTERFACE AND ALL CELLS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to set-recharge-all-patches
  ask patches with [interior-node? = true][set Recharge areal-recharge]
end


to set-recharge-single-patch
  if mouse-down?     ;; reports true or false to indicate whether mouse button is down
    [
      ask patch mouse-xcor mouse-ycor
      [ 
        set Recharge areal-recharge
        set plabel-color black
        set plabel Recharge
      ]
    ]
end


to clear-recharge
  ask patches with [interior-node? = true][set Recharge 0]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; MULTIPLIER FUNCTIONS / RECHARGE AND WELLS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calculate-multiplier-sine-recharge                                            ;; assumes that recharge occurs in (southern hemisphere) winter
    let B-recharge (2 * pi / 365)
    let C-recharge 81.75
    set sine-recharge sin (B-recharge * (ticks - C-recharge) * 180 / pi)         ;; note that the sine function here is in degrees, not in radians
    if sine-recharge < 0 [set sine-recharge 0] 
end


to calculate-multiplier-sine-wells                                               ;; assumes that pumping occurs in (southern hemisphere) summer
    let B-well (2 * pi / 365)
    let C-well -81.75
    set sine-well sin (B-well * (ticks - C-well) * 180 / pi)                     ;; note that the sine function here is in degrees, not in radians
    if sine-well < 0 [set sine-well 0] 
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; SETUP PROCEDURES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-model
  ca
  reset-ticks
  setup-world
  setup-bc
  setup-initial-heads
  setup-hydraulic-parameters
  setup-sources
  setup-maxmin-heads
  setup-neighbours
  setup-patch-numbering
  setup-patches-adjacent-to-boundary
  refresh-view
end


to setup-world
  resize-world 0 (N + 1) 0 (M + 1)                            ;; defines number of cells in X (N) and Y (M) direction, leaves additional patches for the boundary conditions NOTE: N=M
  ask patches [set pcolor white]                              ;; set initial colour of patches to white
  ask patches [set well? false set recharge? false ]           ;; initially none of the patches have wells or recharge, setup-pumping and setup-recharge modify this tag
end


to setup-bc
  ask patches [set no-flow? false set fixed-head? false set interior-node? true]                                                            ;; initialize these location indicators: all nodes are interior and no bc tags - easier for further coding
  
  ask patches [if (pxcor = min-pxcor) or (pxcor = max-pxcor) or (pycor = max-pycor) or (pycor = min-pycor) [set interior-node? false]]      ;; tag interior and boundary nodes
  
  ifelse left-bc = "no-flow" 
  [ask patches with [pxcor = 0][set pcolor black set no-flow? true set interior-node? false set fixed-head? false set K-patch 0 set T-patch 0]]
  [ask patches with [pxcor = 0][set pcolor blue set fixed-head? true set interior-node? false set H left-bc-head]]
  
  ifelse right-bc = "no-flow" 
  [ask patches with [pxcor = max-pxcor][set pcolor black set no-flow? true set interior-node? false set fixed-head? false set K-patch 0 set T-patch 0]]
  [ask patches with [pxcor = max-pxcor][set pcolor blue set fixed-head? true set interior-node? false set H right-bc-head]]
  
  ifelse top-bc = "no-flow" 
  [ask patches with [pycor = max-pycor][set pcolor black set no-flow? true set interior-node? false set fixed-head? false set K-patch 0 set T-patch 0]]
  [ask patches with [pycor = max-pycor][set pcolor blue set fixed-head? true set interior-node? false set H top-bc-head]]
  
  ifelse bottom-bc = "no-flow" 
  [ask patches with [pycor = min-pycor][set pcolor black set no-flow? true set interior-node? false set fixed-head? false set K-patch 0 set T-patch 0]]
  [ask patches with [pycor = min-pycor][set pcolor blue set fixed-head? true set interior-node? false set H bottom-bc-head]]
end


to setup-initial-heads
  ask patches with [interior-node? = true] [set H initial-heads]                               ;; set the initial heads in interior nodes
end


to setup-hydraulic-parameters                                       ;; sets a homogeneous conductivity to the whole model using the input box on the interface (avoids having cells with no K value)
  ask patches with [no-flow? = false]
  [
    set K-patch K 
    ifelse aquifer-type = "confined"
    [set T-patch (K-patch * aquifer-thickness) set S-patch S]       ;; sets transmissivity and storage for confined conditions
    [set T-patch (K-patch * H) set S-patch Sy]                      ;; sets transmissivity and storage for unconfined conditions                    
  ]
end


to setup-sources
  ask patches [set well? false]         ;; initialize patch tags without wells
  ask patches [set injection? false]    ;; initialize patch tags without injection
  ask patches [set recharge? false]     ;; initialize patch tags without recharge
  ask patches [set ET? false]           ;; initialize patch tags without ET
  ask patches [set DRAIN? false]        ;; initialize patch tags without DRAIN
  ask patches [set RIV? false]          ;; initialize patch tags without leakage from river
  
  ;; Sources and stresses can be added directly through the interface, otherwise manually (directly in the code) here:
  
  ;; pumping wells, units of [m3/day] ~ [L3/T]
  ;; injection wells, units of [m3/day] ~ [L3/T]
  ;; recharge cells, units of [L3/L2*T] ~ [L/T]
  ;; ET cells, units of [L3/L2*T] ~ [L/T]
  ;; DRAIN cells, units of [L3/T] 
  ;; RIV cells, units of [L3/T]
                                                      
end


to setup-neighbours
  ask patches with [interior-node? = true]
  [
    set N-neighbour patch-at 0 1                                          ;; neighbour North of this patch
    set S-neighbour patch-at 0 -1                                         ;; neighbour South of this patch
    set E-neighbour patch-at 1 0                                          ;; neighbour East of this patch
    set W-neighbour patch-at -1 0                                         ;; neighbour West of this patch
  ]
end


to setup-patch-numbering               
  ask patches with [interior-node? = true]
  [set patch-number-tag pxcor + (M - pycor) * N]
end


to setup-patches-adjacent-to-boundary
  ask patches with [interior-node? = true]
  [
    ifelse (pxcor = (min-pxcor + 1) or pxcor = (max-pxcor - 1) or pycor = (min-pycor + 1) or pycor = (max-pycor - 1))
    [set adjacent-to-boundary? true]
    [set adjacent-to-boundary? false] 
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; MAIN ITERATION LOOP ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to GO 
  ;; *******INSERT SOCIAL MODEL FUNCTIONS IN THE FOLLOWING BLOCK******
  if social-model = TRUE                         ;; option to deactivate the social model and run groundwater model only
  [
    ask one-of patches [SET-WATER-PRICE]
    ;    ask patches 
    ;    [set-pumping-cost]
    UPDATE-CANAL-WATER-PRICE
    
    if drilling? = TRUE
    [
      if ticks > 0 [UPDATE-MAX-WATER-DEPTH]      ;; keeps track of the max water table depth observed in the basin 
      CHECK-IF-DRY-WELLS                         ;; if farmer well goes dry, farmer hires driller. Drills new well at a cost defined by 'drilling-cost' slider in the interface
      DRILL-A-WELL                               ;; when a farmer contacts a drilling company 
    ]    
    
    ESTABLISH-TREES
    CLAIM-WATER
    ; UPDATE-PUMPING-RATES
  ]
  
  UPDATE-INFLOW-TO-BASIN
  iterate
  ifelse water-level-labels? = FALSE [set view "crops"][set view "heads (values)"]
  refresh-view   
  
  if social-model = TRUE 
    [
      ask turtles[
        set-pumping-costs-ex-post
        CALCULATE-MARGIN]   ;; ex post calculation
      SELL-CROPS
      BANKRUPT
      UPDATE-CROP-COUNTS
      UPDATE-CROP-VIEW
      
      update-crop-distribution
      update-lorenz-and-gini
    ]
  
  tick
  prepare-output
  set year year + 1 
  set plots-ready? TRUE
end


to iterate
  if aquifer-type = "unconfined" [prepare-equations solve-once update-unconfined-transmissivities]
  if aquifer-type = "confined" [solve-once]            
end


to update-unconfined-transmissivities
  ask patches [set T-patch (K-patch * H)]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; POPULATE THE MODEL WITH INITIAL HEADS — STEADY STATE RUN WITH NO STRESSES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to initialize-heads                                                    ;; runs the model ONCE, with all stresses = 0 to set up the initial heads
    reset-ticks

    ask patches
    [
      if aquifer-type = "confined"                                     ;; sets transmissivity and storage for confined conditions
      [set T-patch (K-patch * aquifer-thickness) set S-patch S]       
      
      if aquifer-type = "unconfined"                                   ;; sets transmissivity and storage for unconfined conditions 
      [set T-patch (K-patch * H) set S-patch Sy]                      
      
      set Qwell-temp Qwell                                             ;; save the original values to restore after the init, they are needed for the steady/transient simulation
      set Qinjection-temp Qinjection 
      set Recharge-temp Recharge
      set S-patch-temp S-patch
      set ET-max-rate-temp ET-max-rate-patch
      
      set Qwell 0                                                      ;; we set these parameters to zero for the init, which is steady state without any stresses
      set Qinjection 0 
      set Recharge 0
      set S-patch 0
      set ET-max-rate-patch 0     
    ] 
   
    ask patches with [ET? = true]                                      ;; sets the ET cell to a fixed-head condition to give a reasonalble steady-state solution considering long-term equilibrium at H = ground elevation - extinction depth
    [
      set H (ET-land-surface-elevation - ET-extinction-depth-patch)
      set fixed-head? true
    ]
   
    prepare-equations                                                  ;; set up the equations for this solution
    iterate 
    
    refresh-view 
   
    ;; now prepare equations for the main simulation
    
    if solver = "steady-state" [prepare-equations-steadystate]                                                   ;; note that the fixed-head assumption for ET cells is maintained here
    if solver = "transient" [ask patches with [ET? = true][set fixed-head? false] prepare-equations-transient]   ;; remove the fixed-head assumption of ET cells for subsequent transient simulations
  
    set-initial-heads
end


to prepare-equations-steadystate
  ask patches                         
  [
    set Qwell Qwell-temp 
    set Qinjection Qinjection-temp 
    set Recharge Recharge-temp 
    set S-patch 0                                          ;; S=0 for steady state simulations and stresses are activated
    set ET-max-rate-patch ET-max-rate-temp                       
  ]                      
  prepare-equations
end


to prepare-equations-transient
  ask patches 
  [
    set Qwell Qwell-temp 
    set Qinjection Qinjection-temp 
    set Recharge Recharge-temp 
    set S-patch S-patch-temp 
    set ET-max-rate-patch ET-max-rate-temp
  ]                                                        ;; restores original parameters for the simulation
  prepare-equations
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; PREPARE EQUATIONS AND DEFINE THE INVERSE OF MATRIX A ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; executed ONCE for CONFINED conditions (T does not change during the simulation) and ITERATIVELY for unconfined conditions (T=Kxh varies during the simulation)

to prepare-equations                                                      
  calculate-conductances
  modify-conductances-for-bcs
  build-matrix-A
  remove-inactive-cells-from-A-matrix
  calculate-inverse-A  
end


to calculate-conductances                                                ;; sets the interblock transmissivities using the harmonic mean 
ask patches with [interior-node? = true] 
[
  set TN [T-patch] of N-neighbour                                        ;; T of patch North of here
  set TS [T-patch] of S-neighbour                                        ;; T of patch South of here
  set TE [T-patch] of E-neighbour                                        ;; T of patch East of here
  set TW [T-patch] of W-neighbour                                        ;; T of patch West of here
  set TC [T-patch] of self
 
  ;; using the harmonic mean of conductivity
 
  set Ncof (1 / (delta ^ 2)) * (2 * TC * TN) / (TC + TN)                 ;; Coeff for head North of this patch
  set Scof (1 / (delta ^ 2)) * (2 * TC * TS) / (TC + TS)                 ;; Coeff for head South of this patch
  set Ecof (1 / (delta ^ 2)) * (2 * TC * TE) / (TC + TE)                 ;; Coeff for head East of this patch
  set Wcof (1 / (delta ^ 2)) * (2 * TC * TW) / (TC + TW)                 ;; Coeff for head West of this patch
  set Ccof -1 * (Ncof + Scof + Ecof + Wcof + S-patch / delta-t)          ;; Coeff for head at this patch
]
end


to modify-conductances-for-bcs
ask patches with [interior-node? = true]
[
  ;; conductance to/from no-flow boundaries are zero
  
  if [no-flow?] of N-neighbour = true [set Ncof 0]
  if [no-flow?] of S-neighbour = true [set Scof 0]
  if [no-flow?] of E-neighbour = true [set Ecof 0]
  if [no-flow?] of W-neighbour = true [set Wcof 0]
  
  set Ccof -1 * (Ncof + Scof + Ecof + Wcof + (S-patch / delta-t))                ;; needed to reflect any of the above changes 
]
end


to build-matrix-A                                                                ;; AxB=C ~ NOTE: N=M
  set number-of-unknowns N * M
  
  ;set C matrix:make-constant (number-of-unknowns) 1 0                           ;; a column vector with 0's: here goes the SourceTerm + terms derived from bc's
  set A matrix:make-constant number-of-unknowns number-of-unknowns 0             ;; a N^2xM^2 matrix of 0's: here we will add the equation coeffs
             
  ask patches with [interior-node? = true][set interior-but-corner? false]       ;; set all the interior nodes to false, we tag the corner ones later          
                                                           
  ;; COMPLETE THE "A" MATRIX

  ;; interior patches not adjacent to a boundary have all the the terms
  ask patches with [interior-node? = true and adjacent-to-boundary? = false]
  [
    let my-position ([patch-number-tag] of self - 1)               ;; NOTE: indexing of the matrix extension starts at 0, not 1
    
    matrix:set A (my-position) (my-position) [Ccof] of self        ;; set the center term (diagonal) 
    matrix:set A (my-position) (my-position - N) [Ncof] of self    ;; set the north term
    matrix:set A (my-position) (my-position - 1) [Wcof] of self    ;; set the west term
    matrix:set A (my-position) (my-position + 1) [Ecof] of self    ;; set the east term
    matrix:set A (my-position) (my-position + N) [Scof] of self    ;; set the south term
  ]
  
  ;; top left corner
  ask patches with [interior-node? = true and pycor = max-pycor - 1 and pxcor = min-pxcor + 1]
  [
    let my-position ([patch-number-tag] of self - 1)               ;; NOTE: indexing of the matrix extension starts at 0, not 1 
   
    matrix:set A (my-position) (my-position) [Ccof] of self        ;; set the center term (diagonal) 
    matrix:set A (my-position) (my-position + 1) [Ecof] of self    ;; set the east term
    matrix:set A (my-position) (my-position + N) [Scof] of self    ;; set the south term
    
    set interior-but-corner? true
  ]
  
  ;; top right corner
  ask patches with [interior-node? = true and pycor = max-pycor - 1 and pxcor = max-pxcor - 1]
  [
    let my-position ([patch-number-tag] of self - 1)               ;; NOTE: indexing of the matrix extension starts at 0, not 1 
   
    matrix:set A (my-position) (my-position) [Ccof] of self        ;; set the center term (diagonal) 
    matrix:set A (my-position) (my-position - 1) [Wcof] of self    ;; set the west term
    matrix:set A (my-position) (my-position + N) [Scof] of self    ;; set the south term
    
    set interior-but-corner? true
  ]
  
  ;; bottom left corner
  ask patches with [interior-node? = true and pycor = min-pycor + 1 and pxcor = min-pxcor + 1]
  [
    let my-position ([patch-number-tag] of self - 1)               ;; NOTE: indexing of the matrix extension starts at 0, not 1 
   
    matrix:set A (my-position) (my-position) [Ccof] of self        ;; set the center term (diagonal) 
    matrix:set A (my-position) (my-position - N) [Ncof] of self    ;; set the north term
    matrix:set A (my-position) (my-position + 1) [Ecof] of self    ;; set the east term 
    
    set interior-but-corner? true
  ]
  
  ;; bottom right corner
  ask patches with [interior-node? = true and pycor = min-pycor + 1 and pxcor = max-pxcor - 1]
  [
    let my-position ([patch-number-tag] of self - 1)               ;; NOTE: indexing of the matrix extension starts at 0, not 1 
   
    matrix:set A (my-position) (my-position) [Ccof] of self        ;; set the center term (diagonal) 
    matrix:set A (my-position) (my-position - N) [Ncof] of self    ;; set the north term
    matrix:set A (my-position) (my-position - 1) [Wcof] of self    ;; set the west term    
    
    set interior-but-corner? true
  ]
  
  ;; top row
  ask patches with [interior-node? = true and pycor = max-pycor - 1 and interior-but-corner? = false]
  [
    let my-position ([patch-number-tag] of self - 1)               ;; NOTE: indexing of the matrix extension starts at 0, not 1 
    
    matrix:set A (my-position) (my-position) [Ccof] of self        ;; set the center term (diagonal) 
    matrix:set A (my-position) (my-position - 1) [Wcof] of self    ;; set the west term
    matrix:set A (my-position) (my-position + 1) [Ecof] of self    ;; set the east term
    matrix:set A (my-position) (my-position + N) [Scof] of self    ;; set the south term   
  ]
  
  ;; bottom row
    ask patches with [interior-node? = true and pycor = min-pycor + 1 and interior-but-corner? = false]
  [
    let my-position ([patch-number-tag] of self - 1)               ;; NOTE: indexing of the matrix extension starts at 0, not 1
    
    matrix:set A (my-position) (my-position) [Ccof] of self        ;; set the center term (diagonal) 
    matrix:set A (my-position) (my-position - N) [Ncof] of self    ;; set the north term
    matrix:set A (my-position) (my-position - 1) [Wcof] of self    ;; set the west term
    matrix:set A (my-position) (my-position + 1) [Ecof] of self    ;; set the east term   
  ]
  
  ;; left column
    ask patches with [interior-node? = true and pxcor = min-pxcor + 1 and interior-but-corner? = false]
  [
    let my-position ([patch-number-tag] of self - 1)               ;; NOTE: indexing of the matrix extension starts at 0, not 1
    
    matrix:set A (my-position) (my-position) [Ccof] of self        ;; set the center term (diagonal) 
    matrix:set A (my-position) (my-position - N) [Ncof] of self    ;; set the north term
    matrix:set A (my-position) (my-position + 1) [Ecof] of self    ;; set the east term
    matrix:set A (my-position) (my-position + N) [Scof] of self    ;; set the south term
  ]
  
  ;; right column
    ask patches with [interior-node? = true and pxcor = max-pxcor - 1 and interior-but-corner? = false]
  [
    let my-position ([patch-number-tag] of self - 1)               ;; NOTE: indexing of the matrix extension starts at 0, not 1
    
    matrix:set A (my-position) (my-position) [Ccof] of self        ;; set the center term (diagonal) 
    matrix:set A (my-position) (my-position - N) [Ncof] of self    ;; set the north term
    matrix:set A (my-position) (my-position - 1) [Wcof] of self    ;; set the west term
    matrix:set A (my-position) (my-position + N) [Scof] of self    ;; set the south term 
  ]
end


to calculate-inverse-A
  set inverse-A matrix:inverse A
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; ITERATION SUBROUTINES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to solve-once
  calculate-source-term                       ;; comprises the sources (wells, ET) and sinks (injection, recharge)  
  calculate-RHS                               ;; obtain the RHS of each equation
  build-matrix-C                              ;; these only need to be calculated at the beggining of each iteration
  reduce-matrix-C                             ;; eliminate the sourceterms corresponding to no-flow and fixed-head cells
  solve-system-of-equations                   ;; solve the system of equations AxB=C
  list-of-active-cells                        ;; lists the cells to be updated with new heads
  extract-Hinew                               ;; the new heads are mapped into the Hinew variable of the grid
  update-heads-to-grid                        ;; H <-- Hinew / the iteration has ended so we update the H values for the iterations corresponding to the next timestep
end


to calculate-source-term
  
  ifelse (sine-recharge-multiplier? = true and solver = "transient") [calculate-multiplier-sine-recharge][set sine-recharge 1]
  ifelse (sine-well-multiplier? = true and solver = "transient") [calculate-multiplier-sine-wells][set sine-well 1]
  
  if solver = "steady-state"                                      ;; IMPORTANT NOTE: when using ET, DRAIN and RIV (head-dependent) cells, set initial conditions using a long-term transient simulation
  [
  ]                
  
  if solver = "transient"                                         ;; recalculate the ET and DRAIN terms before each iteration   
  [
    ask patches with [ET? = true][calculate-ET-discharge]
    ask patches with [DRAIN? = true][calculate-DRAIN-discharge]
    ask patches with [RIV? = true][calculate-RIV-term]
  ]                                                                                             
  
  ;;FOR THE CELLS THAT ARE NOT FARMER WELLS I.E. TOWN WELLS PUMP CONSTANTLY THROUGHOUT THE YEAR
  ask patches with [interior-node? = true]
  [
    set SourceTerm ((Qwell - Qinjection) / (delta ^ 2)) + (- Recharge * sine-recharge) + (ET) + (DRAIN / (delta ^ 2)) + (RIV / (delta ^ 2))      ;; units consistent [L3/L2*T] ~ [L/T]
    
    ;; look at neighbours for fixed-head boundary, we know the heads here and therefore move this term to the RHS
  
    if [fixed-head?] of N-neighbour = true [set SourceTerm (SourceTerm - ([H] of N-neighbour) * (Ncof))]
    if [fixed-head?] of S-neighbour = true [set SourceTerm (SourceTerm - ([H] of S-neighbour) * (Scof))]
    if [fixed-head?] of E-neighbour = true [set SourceTerm (SourceTerm - ([H] of E-neighbour) * (Ecof))]
    if [fixed-head?] of W-neighbour = true [set SourceTerm (SourceTerm - ([H] of W-neighbour) * (Wcof))] 
  ]
  
  ;;MAKES FARMERS PUMP USING A SEASONAL PATTERN
  ask patches with [(interior-node? = true and farmer-well-here? = true) or (interior-node? = true and well? = true)]
  [
    set SourceTerm ((Qwell - Qinjection) * sine-well / (delta ^ 2)) + (- Recharge * sine-recharge) + (ET) + (DRAIN / (delta ^ 2)) + (RIV / (delta ^ 2))      ;; units consistent [L3/L2*T] ~ [L/T]
    
    ;; look at neighbours for fixed-head boundary, we know the heads here and therefore move this term to the RHS
  
    if [fixed-head?] of N-neighbour = true [set SourceTerm (SourceTerm - ([H] of N-neighbour) * (Ncof))]
    if [fixed-head?] of S-neighbour = true [set SourceTerm (SourceTerm - ([H] of S-neighbour) * (Scof))]
    if [fixed-head?] of E-neighbour = true [set SourceTerm (SourceTerm - ([H] of E-neighbour) * (Ecof))]
    if [fixed-head?] of W-neighbour = true [set SourceTerm (SourceTerm - ([H] of W-neighbour) * (Wcof))] 
  ]
  
  ;;FIXED FLUX CELLS 
  ask patches with [interior-node? = true and injection? = true and fixed-flux? = true]
  [
    set SourceTerm ((Qwell - Qinjection) / (delta ^ 2)) + (- Recharge * sine-recharge) + (ET) + (DRAIN / (delta ^ 2)) + (RIV / (delta ^ 2))      ;; units consistent [L3/L2*T] ~ [L/T]
    
    ;; look at neighbours for fixed-head boundary, we know the heads here and therefore move this term to the RHS
  
    if [fixed-head?] of N-neighbour = true [set SourceTerm (SourceTerm - ([H] of N-neighbour) * (Ncof))]
    if [fixed-head?] of S-neighbour = true [set SourceTerm (SourceTerm - ([H] of S-neighbour) * (Scof))]
    if [fixed-head?] of E-neighbour = true [set SourceTerm (SourceTerm - ([H] of E-neighbour) * (Ecof))]
    if [fixed-head?] of W-neighbour = true [set SourceTerm (SourceTerm - ([H] of W-neighbour) * (Wcof))] 
  ] 
end


to calculate-RHS
  ask patches with [interior-node? = true]
  [set RHS (-1 * S-patch) / delta-t * (H) + SourceTerm]
end


to build-matrix-C
  set C matrix:make-constant (number-of-unknowns) 1 0  
  ask patches with [interior-node? = true]
  [
    let my-position ([patch-number-tag] of self - 1)                ;; NOTE: indexing of the matrix extension starts at 0, not 1
    matrix:set C (my-position) 0 [RHS] of self                      ;; set Q term of the RHS of the equation   
  ]
end


to solve-system-of-equations
  set solution-vector matrix:times inverse-A C
end


to extract-Hinew
  let i 0
  while [i <= length active-cells-remap - 1]
  [
    ask patches with [patch-number-tag = (item i active-cells-remap)][set Hinew matrix:get (solution-vector) i 0]
    set i i + 1
  ]
end


to update-heads-to-grid
  ask patches with [interior-node? = true and no-flow? = false and fixed-head? = false][set H Hinew]                    ;; replace the old head value with the new one
end


to calculate-ET-discharge                                                                                             ;; recalculates the ET term at each iteration
  ask patches with [ET? = true]
  [
   if H > ET-land-surface-elevation-patch [set ET ET-max-rate-patch]                                                  ;; if head is greater than land surface elevation
   
   if H >= (ET-land-surface-elevation-patch - ET-extinction-depth-patch)                                              ;; if head is between land surface and extinction depth
   and H <= (ET-land-surface-elevation-patch) 
   [set ET ET-max-rate-patch * (H - (ET-land-surface-elevation-patch - ET-extinction-depth-patch)) / ET-extinction-depth-patch]   ;; a linear relationship between the max rate and zero
   
   if H < (ET-land-surface-elevation-patch - ET-extinction-depth-patch) [set ET 0]                                    ;; if head is lower than extinction depth     
  ]
  
end


to calculate-DRAIN-discharge                                                                                          ;; recalculates the springflow at each iteration NOTE: this values is in m3/day, needs to be divided by the cell area in the sourceterm equation
  ask patches with [DRAIN? = true]
  [
   if H > DRAIN-elevation-patch [set DRAIN DRAIN-conductance-patch * (H - DRAIN-elevation-patch)]
   if H <= DRAIN-elevation-patch [set DRAIN 0]
  ]
end

to calculate-RIV-term                                                                                         ;; recalculates the springflow at each iteration NOTE: this values is in m3/day, needs to be divided by the cell area in the sourceterm equation
  ask patches with [RIV? = true]
  [
   if H > RIV-elevation-patch [set RIV RIV-conductance-patch * (H - RIV-elevation-patch)]
   if H <= RIV-elevation-patch and H > RIV-bottom-patch [set RIV RIV-conductance-patch * (RIV-elevation-patch - H) * -1]
   if H <= RIV-elevation-patch and H < RIV-bottom-patch [set RIV RIV-conductance-patch * (RIV-elevation-patch - RIV-bottom-patch) * -1]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; REMOVE INACTIVE CELLS FROM CALCULATIONS / SPEEDS UP SOLUTION ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to remove-inactive-cells-from-A-matrix
  list-inactive-cells
  define-rows-to-remove
  remove-rows
  transform-A-from-row-list-to-matrix
  remove-columns
  transform-A-from-column-list-to-matrix
  set A A-reduced
end


to remove-rows
  transform-A-to-row-list 
  remove-rows-from-A 
end


to remove-columns
  transform-A-to-column-list
  remove-columns-from-A
end


to reduce-matrix-C
  transform-C-to-row-list
  remove-rows-from-C
  transform-C-from-row-list-to-matrix
  set C C-reduced
end


to list-inactive-cells                                                              ;; these are either no-flow or fixed head cells in the model. We need to remove these from the calculations
  ask patches with [interior-node? = true]
  [
    set no-flow-cells-list (list [patch-number-tag] of patches with [no-flow? = true and interior-node? = true])
    set fixed-head-cells-list (list [patch-number-tag] of patches with [fixed-head? = true and interior-node? = true])
    set remove-cells-list sentence no-flow-cells-list fixed-head-cells-list
    set remove-cells-list sort sentence item 0 remove-cells-list item 1 remove-cells-list
  ]
end


to define-rows-to-remove                                                            ;; creates a list containing the index of patches that we will remove. Note that the indexing begints at 0. It considers the fact that the remove-item function removes one row at a time and the indexes change during this process
  let i 0
  set remove-cells-list-corrected []
  while [i <= (length remove-cells-list - 1) ]
  [
    let index-of-cell-to-remove (item i remove-cells-list) - i - 1
    set remove-cells-list-corrected lput index-of-cell-to-remove remove-cells-list-corrected
    set i i + 1
  ]
end


to transform-A-to-row-list                                                         ;; transforms matrix A into list form (rows)
  set A-row-list matrix:to-row-list A
end


to remove-rows-from-A                                                              ;; removes the rows one by one                      
  let i 0
  set A-row-reduced-list A-row-list
  while [i <= (length remove-cells-list - 1) ]
  [
    let index-of-cell-to-remove (item i remove-cells-list-corrected)
    set A-row-reduced-list remove-item index-of-cell-to-remove A-row-reduced-list
    set i i + 1
  ]
end

to transform-A-from-row-list-to-matrix                                              ;; back to matrix form
  set A-row-reduced-matrix matrix:from-row-list A-row-reduced-list
end


to transform-A-to-column-list                                                      ;; transforms matrix A into list form (columns)
  set A-column-list matrix:to-column-list A-row-reduced-matrix
end


to remove-columns-from-A                                                           ;; removes the rows one by one                      
  let i 0
  set A-column-reduced-list A-column-list
  while [i <= (length remove-cells-list - 1) ]
  [
    let index-of-cell-to-remove (item i remove-cells-list-corrected)
    set A-column-reduced-list remove-item index-of-cell-to-remove A-column-reduced-list
    set i i + 1
  ]
end


to transform-A-from-column-list-to-matrix                                              ;; back to matrix form
  set A-reduced matrix:from-column-list A-column-reduced-list
end


to transform-C-to-row-list
   set C-row-list matrix:to-row-list C 
end
    
    
to remove-rows-from-C
  let i 0
  set C-row-reduced-list C-row-list
  while [i <= (length remove-cells-list - 1) ]
  [
    let index-of-cell-to-remove (item i remove-cells-list-corrected)
    set C-row-reduced-list remove-item index-of-cell-to-remove C-row-reduced-list
    set i i + 1
  ]  
end
  
  
to transform-C-from-row-list-to-matrix
    set C-reduced matrix:from-row-list C-row-reduced-list    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; ACTIVE CELLS FOR MAPPING SOLUTION VECTOR ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to list-of-active-cells
  set active-cells-remap (list [patch-number-tag] of patches with [no-flow? = false and fixed-head? = false and interior-node? = true])
  set active-cells-remap sort item 0 active-cells-remap
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; WRITE HEADS TO FILE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to write-output-heads
  file-delete "heads.txt"
  file-open "heads.txt"
  let j 1
  while [j <= N]
  [
    let i 1
    while [i <= M]
    [
      if i <= M [ask patch i j [file-write H]]
      if i = M [ask patch i j [file-print ""]]
      set i i + 1
    ]  
    set j j + 1 
  ]
  file-close
  file-flush
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; BUDGET CALCULATIONS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;WATER BUDGET FOR RIVER
to calculate-river-budget                   ;;call this procedure within the 'go' process to calculate this budget term
  ask patches with [RIV? = true]
  [
    set riv-budget riv-budget + RIV      
  ] 
end


;;TOTAL WATER PUMPED
to calculate-well-budget                   ;;call this procedure within the 'go' process to calculate this budget term
    ask patches with [well? = true]
  [
    set well-budget well-budget + Qwell * sine-well
  ] 
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; RECYCLE CODE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; SOCIAL MODEL PROCEDURES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; SET VALUES & LISTS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; units for 100 hectares... 247 acres

to DEFINE-LISTS
  ask patches [
    set gmlist[]
    set gmlist-temp[]
    set past-crop-index-list []
    set past-crop-list []
    set past-gm-list []
    set past-water-price-list []
  ]  
end 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; MAIN SETUP ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to SETUP-SOCIAL
  ;set-prices
  ;set-yield
  ;set-water-demand
  ;set-other-costs
  ;prepare-GIS-map
  UPDATE-CANAL-WATER-PRICE
  SET-INITIAL-WATER-PRICE
  DEFINE-LISTS
  read-crop-data
  set max-depth-list []
  SETUP-FARMERS
  set model-ready? true
  update-lorenz-and-gini
  reset-ticks
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; SETUPS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to SETUP-FARMERS
  
  ask patches with [(pxcor < min-pxcor + buffer-zone) or (pxcor > max-pxcor - buffer-zone) or (pycor < min-pycor + buffer-zone) or (pycor > max-pycor - buffer-zone)] [set farmer-here? false]
  ask patches with [(pxcor >= min-pxcor + buffer-zone) and (pxcor <= max-pxcor - buffer-zone) and (pycor >= min-pycor + buffer-zone) and (pycor <= max-pycor - buffer-zone)] [set farmer-here? true]
  ask patches with [gw-bank-here? = true] [set farmer-here? false]
  
  ask patches with [farmer-here? = true] 
  [
    set Hground 500                                                    ;; [meters]
    ;set Hinitial 400                                                   ;; [meters]
    set crop-assigned? false
    set bankrupt? false
    set wait-counter 0
    set canal-rights-rank (count patches with [canal-rights-rank > 0] + 1)
    ; alters the distribution of starting-wealth
    
    ifelse starting-wealth-distribution-uniform? 
    [set money random 2 * median-starting-farmer-money]
    ; standard-deviation-size is a proportion of starting-farmer-money
    [set money random-normal median-starting-farmer-money standard-deviation-size * median-starting-farmer-money]
    
    set dry-well? false 
    set tree-grower? false
    set new-well-needed? false
    set well-depth initial-well-depth              ;; sets initial well depths at using the interface slider
    
                                                   ;; what is the farmers decisions making heuristic?
    
    SET-FARMING-STRATEGY
  ]

     
  ;; AEWSD crops: Grains  Alfalfa  Cotton  Onion  Truck  Deciduous  Grapes  Potatoes  Melon  Carrots  Citrus — (SCHUCK ET AL. 2002)
  
  set %-of-Grains item 0 acreage-list-pct
  set %-of-Alfalfa item 1 acreage-list-pct
  set %-of-Cotton item 2 acreage-list-pct
  set %-of-Onion item 3 acreage-list-pct
  set %-of-Truck item 4 acreage-list-pct
  set %-of-Deciduous item 5 acreage-list-pct
  set %-of-Grapes item 6 acreage-list-pct
  set %-of-Potatoes item 7 acreage-list-pct
  set %-of-Melon item 8 acreage-list-pct
  set %-of-Carrots item 9 acreage-list-pct
  set %-of-Citrus item 10 acreage-list-pct
  set %-of-Fallow 0
  
  ask n-of (((N * M) - count patches with [farmer-here? = false]) * %-of-Grains) patches with [farmer-here? = true and crop-assigned? = false]
    [
      set current-crop "Grains"
      set pcolor 32
      set crop-assigned? true
    ]
  
  ask n-of (((N * M) - count patches with [farmer-here? = false]) * %-of-Alfalfa) patches with [farmer-here? = true and crop-assigned? = false]
    [
      set current-crop "Alfalfa"
      set pcolor 64
      set crop-assigned? true
    ]
  
  ask n-of (((N * M) - count patches with [farmer-here? = false]) * %-of-Cotton) patches with [farmer-here? = true and crop-assigned? = false]
    [
      set current-crop "Cotton"
      set pcolor 39
      set crop-assigned? true
    ]
  
  ask n-of (((N * M) - count patches with [farmer-here? = false]) * %-of-Onion) patches with [farmer-here? = true and crop-assigned? = false]
    [
      set current-crop "Onion"
      set pcolor 117
      set crop-assigned? true
    ]
  
  ask n-of (((N * M) - count patches with [farmer-here? = false]) * %-of-Truck) patches with [farmer-here? = true and crop-assigned? = false]
    [
      set current-crop "Truck"
      set pcolor 56
      set crop-assigned? true
    ]
  
  ask n-of (((N * M) - count patches with [farmer-here? = false]) * %-of-Deciduous) patches with [farmer-here? = true and crop-assigned? = false]
    [
      set current-crop "Deciduous"
      set pcolor 68
      set crop-assigned? true
    ]
  
  ask n-of (((N * M) - count patches with [farmer-here? = false]) * %-of-Grapes) patches with [farmer-here? = true and crop-assigned? = false]
    [
      set current-crop "Grapes"
      set pcolor 125
      set crop-assigned? true
    ]
      
  ask n-of (((N * M) - count patches with [farmer-here? = false]) * %-of-Potatoes) patches with [farmer-here? = true and crop-assigned? = false]
    [
      set current-crop "Potatoes"
      set pcolor 37
      set crop-assigned? true
    ]
  
  ask n-of (((N * M) - count patches with [farmer-here? = false]) * %-of-Melon) patches with [farmer-here? = true and crop-assigned? = false]
    [
      set current-crop "Melon"
      set pcolor 15
      set crop-assigned? true
    ]
  
  ask n-of (((N * M) - count patches with [farmer-here? = false]) * %-of-Carrots) patches with [farmer-here? = true and crop-assigned? = false]
    [
      set current-crop "Carrots"
      set pcolor 25
      set crop-assigned? true
    ]
  
  ask n-of (((N * M) - count patches with [farmer-here? = false]) * %-of-Citrus) patches with [farmer-here? = true and crop-assigned? = false]
    [
      set current-crop "Citrus"
      set pcolor 27
      set crop-assigned? true
    ]
  
  UPDATE-CROP-VIEW
end

;;; ******I AM NOT TOO SURE THIS THE THE CORRECT WAY TO ASSIGN AGENTS TO STRATEGIES...*****
;;; % COPYCATS AND % OPTIMIZERS SHOULD BE RELATED, THERE IS NO RELATION IN THE MATH BELOW


to SET-FARMING-STRATEGY
     ifelse random-float 1 < %-optimizers 
      [ set optimizer? true
        set copycat? false
      ] ;;;;
      [ifelse random-float 1 < remaining-%-copycats[
        set optimizer? false
        set copycat? true]
      [set optimizer? false
        set copycat? false]
    ]
  
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; SET-WATER-PRICE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to SET-INITIAL-WATER-PRICE
  ;; ground-water starts out at the same price as canal water
  set water-price canal-water-price
  if canal-water-price = 0 [set water-price 500]
end

to SET-WATER-PRICE
  ifelse prices = true [
    
    ;; mean initial is average level for max-capacity of basin    
    let mean-initial mean [Hground] of patches with [interior-node? = true] 
    
    ;; mean-current is current average level of water in basin 
    let mean-current mean [H] of patches with [interior-node? = true ] 
    
    ;; if ground water is more expensive or the same price as canal water...
    
    
    ;; if using the algorithm does not produce a negative price of ground-water  
    if water-price  + price-sensitivity * ln (target-water-level * mean-initial / mean-current) >= 0 [
;      ifelse (count patches with [current-crop = "fallow" and interior-node? = true and injection? = false]) / (count patches with [interior-node? = true and injection? = false]) < .4 
      
       
        
        ;; if level is 5 percent above target, return ground water price to canal price - 1
        ;; this will cause farmers to substitute toward ground water until the total cost of use if equal to or less than canal water.
        ;; ensures that ground water does not remain overpriced
        ifelse water-price >= canal-water-price or (target-water-level * mean-initial / mean-current) > 0[
          set water-price water-price + price-sensitivity * ln (target-water-level * mean-initial / mean-current)
          if (mean-current >= 1.05 * target-water-level * mean-initial and water-price  + price-sensitivity * ln (target-water-level * mean-initial / mean-current) >= 0)
           or mean-current <= 95. * target-water-level * mean-initial 
          [set water-price  water-price +  price-sensitivity  * ln (target-water-level * mean-initial / mean-current)
;            if (mean-current >= 1.05 * target-water-level * mean-initial and water-price  + price-sensitivity * ln (target-water-level * mean-initial / mean-current) >= 0)
;            or mean-current <= .95 * target-water-level * mean-initial 
;            [set water-price  water-price +  price-sensitivity  * ln (target-water-level * mean-initial / mean-current)]
          ]      
          
          ;; otherwise set the current price according to the basic algorithm
          
          ;/ (target-water-level * mean-initial)
          
        ]
        
        
        
        [ set water-price water-price + price-sensitivity / 4 * ln (target-water-level * mean-initial / mean-current)]
      ]

      ;[]
      
    ]
  ;]
  
  
  
  [set water-price 0]
  
  
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; FARMER PROCEDURES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; if trees, farmer must drill to avoid loss, if unable to drill any deeper, then go fallow or switch crop
to ESTABLISH-TREES
  ask patches with [farmer-here? = true and tree-grower? = true]
  [
    set wait-counter wait-counter + 1    
    if wait-counter > item current-crop-index waiting-period-list [set tree-grower? false]  ;; no longer waiting for tree yield                    
  ]
end

to CHOOSE-CROPS
  set mutate? false
  if  random-float 1 < mutation-rate [set mutate? true]
  if optimizer? [optimize]
  if copycat? [copy]
  if copycat? = false and optimizer? = false and mutate? = true                    ;;; ******* THESE CONDITIONS MIGHT CONTRADICT EACH OTHER ... ***********
  [
    set current-crop one-of crop-list
    set current-crop-index position current-crop crop-list  
  ]
  if current-crop-temp != current-crop [
    set wait-counter 0
  ]
  ifelse current-crop-index < 5 and wait-counter < item current-crop-index waiting-period-list [
    set tree-grower? true
  ]
  [
    set tree-grower? false
  ]
  if interior-node? = true and dry-well? = true and canal-water? = false[ 
    set current-crop "fallow"
    set current-crop-index position current-crop crop-list
  ]
end

to CALCULATE-MARGIN 

  ;  FOR REFERENCE 
  ;  set cropnames-list item 0 fileList
  ;  set acreage-list item 1 fileList
  ;  set yield-list item 2 fileList
  ;  set wdAF-list item 3 fileList
  ;  set wdM3-list item 4 fileList
  ;  set netprice-list item 5 fileList
  
  

  ;; if canal water is cheaper than pumping well water or if no well-water available
  ;; and if canal water is available: use canal water
  ifelse ((canal-water-price <= (pumping-cost)) or (dry-well? = true))  and canal-water-used < max-canal-flow * drought-index and max-canal-flow > 0 [
     
    set canal-water? true
    
    ;; Use canal water, trees are established
      (foreach netprice-list wdM3-list[ ; oc-list [
          set gmlist-temp lput 
        (?1  - (canal-water-price * (?2 / (1233.48)) )) gmlist-temp]) ;; "" / 1233.48 converts to acre-feet
      set gmlist gmlist-temp
      set gmlist-temp[]
      
   
  ]
  
  ;;pumping cost includes price of groundwater. If prices are off, Pgw = 0
  
  [
   ;; Use groundwater
    
   set canal-water? false
   if dry-well? = false[
    (foreach netprice-list wdM3-list [set gmlist-temp lput 
      (?1  - ((pumping-cost)* ?2 / (1233.48) )) gmlist-temp]) ;; pumping-cost + water-price = total cost of ground water use per unit, converted to acre-feet
      set gmlist gmlist-temp
      set gmlist-temp[]
    ]
    ]
end

to OPTIMIZE
  set current-crop-temp current-crop      
  ifelse bankrupt? = false and item (max-item gmlist) gmlist > 0 [
    
    ;; if not bankrupt, choose crop that yields highest return. Assume optimizers do "research"
    set current-crop (item (max-item gmlist) crop-list)
    set current-crop-index position current-crop crop-list  
  ]
  [
    
    ;; if bankrupt or gmlist elements are all negative, do not plant
    set current-crop "fallow"
    set current-crop-index position current-crop crop-list
    set tree-grower? false
  ]
end



to COPY
  set current-crop-temp current-crop
  let rich-neighbor max-one-of neighbors with [farmer-here? = true][money]
  ifelse bankrupt? = false [
    
    ifelse (mutate? = false and [money] of rich-neighbor >= 0 )
      [
        ;; If not bankrupt, copy your richest neighbor's crop
        set current-crop [current-crop] of rich-neighbor
        set current-crop-index [current-crop-index] of rich-neighbor
      ]
      [ ;; At some rate - mutation-rate - plant a radnom crop. We will improve agent intlligence by choosing one of top-yielding crops
        
        set current-crop-index random 12                    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; choose from top 5 crops
        set current-crop item current-crop-index crop-list  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ]                                                     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ]
  ;; if bankrupt, do not grow
    [
      set current-crop "fallow"
      set current-crop-index position current-crop crop-list
      set tree-grower? false
    ]
end


to BANKRUPT
  ask patches with [farmer-here? = true] [;    [injection? = false and interior-node? = true and tree-grower? = false ][
    set bankrupt? false
    if money < 0 [
      set money 0
      set current-crop "fallow"
      set current-crop-index position current-crop crop-list
      set bankrupt? true
      set wait-counter 0
      if change-strategy-if-bankrupt? = true[SET-FARMING-STRATEGY]
      
    ]
  ]
  ;; global for counting bankruptcies
  set count-bankrupt count-bankrupt + (count patches with [bankrupt? = true])
end

to RECORD-HISTORY
  
  set past-crop-index-list fput current-crop-index past-crop-index-list
  set past-crop-list fput item current-crop-index crop-list past-crop-list
  set past-gm-list fput item current-crop-index gmlist past-gm-list
  ifelse canal-water? = true 
  [set past-water-price-list fput canal-water-price past-water-price-list]
  [set past-water-price-list fput pumping-cost past-water-price-list]
  
  if length past-crop-index-list > 10 [
    set past-crop-index-list remove-item 10 past-crop-index-list
    set past-crop-list remove-item 10 past-crop-list
    set past-gm-list remove-item 10 past-gm-list
    set past-water-price-list remove-item 10 past-water-price-list
  ]
end

to CHECK-PAST-REVENUES
  ; Adaptive epectations turn on if farmer earned losses for past 3 periods
  let max-past-gm-index max-item past-gm-list
  if item 0 past-gm-list < 0[
    if item 1 past-gm-list < 0[
      ;     if item 2 past-gm-list < 0[
      
      ;; if max revenue was greater than 0, grow crop that generate that return
      ifelse item (max-past-gm-index) past-gm-list > 0 and item (max-past-gm-index) past-crop-list != item 0 past-crop-list and item (max-past-gm-index) past-crop-list != item 1 past-crop-list 
      [; and max-past-gm-index  [;and max-past-gm-index != 1 and max-past-gm-index != 2[
        ifelse canal-water? = true[ 
          ;  ifelse item (max-past-gm-index) past-water-price-list >= canal-water-price [
          set current-crop item (max-past-gm-index) past-crop-list
          set current-crop-index position current-crop crop-list
          ;set gmlist item current-crop-index gmlist
          ;  ]
          
          ;  [
          ;          set current-crop "fallow"
          ;          set current-crop-index position current-crop crop-list
          ;set gmlist item current-crop-index gmlist
          ;    ]
        ]
        [ 
          ; ifelse item (max-past-gm-index) past-water-price-list >= pumping-cost [
          set current-crop item (max-past-gm-index) past-crop-list
          set current-crop-index position current-crop crop-list
          ;set gmlist item current-crop-index gmlist
          ; ]
          
          ;[
          ;          set current-crop "fallow"
          ;          set current-crop-index position current-crop crop-list
          ;set gmlist item current-crop-index gmlist
          ;  ]
        ]
        
      ]
      [ set current-crop "fallow"
        set current-crop-index position current-crop crop-list]
    ]];]
end

to CLAIM-WATER
  ;; MAIN SCRIPT, calls functions for eac
  ;; agents run serially according to canal-rights-rank
  
  let i 0
  
  ;; Globals that keep track of total canal and groundwater used in each period
  set canal-water-used 0
  set well-water-used 0
  
  
  ;; Agent action executed according to canal-rights-rank
  while [i <= count patches with [farmer-here? = true]]
  [ask patches with [canal-rights-rank = i and farmer-here? = true] [
    ;; set water use to 0 at beginning of period    
    set Qcanal 0
    set Qwell 0
    CHECK-IF-DRY-WELLS
    if drilling? and money > drilling-cost[DRILL-A-WELL]
    REPAY-WELL
    set-pumping-cost-ex-ante ;; make expected pumping cost, farmer can only see what the level currently is, not what the level will be
    ;; Here farmer will choose between canal and groundwater, also identifies whether trees are being established or are mature
    CALCULATE-MARGIN
    
    ;; crops chosen either by optimization, copying, or randomly
    CHOOSE-CROPS
    if length past-crop-index-list > 3 and mutate? = false [CHECK-PAST-REVENUES]
    RECORD-HISTORY
      
    ;; under calculate-margin, boolean canal-water? is marked true or false
    ifelse canal-water?
    [UPDATE-CANAL-WATER]
    [UPDATE-PUMPING-RATES]
  ]    
  set i i + 1
  
  ]
end

to UPDATE-CANAL-WATER
  
  ;; No delta-t for canal water, rate of flow is yearly, not daily
 
        if current-crop = "Grains" 
        [set Qcanal item 0 wdM3-list  ]   
        
        if current-crop = "Alfalfa"
        [set Qcanal item 1 wdM3-list  ]
        
        if current-crop = "Cotton" 
        [set Qcanal item 2 wdM3-list  ]
        
        if current-crop = "Onion" 
        [set Qcanal item 3 wdM3-list  ]
        
        if current-crop = "Truck" 
        [set Qcanal item 4 wdM3-list  ]
        
        
        if current-crop = "Deciduous" 
        [set Qcanal item 5 wdM3-list  ]
        
        if current-crop = "Grapes" 
        [set Qcanal item 6 wdM3-list  ]
        
        if current-crop = "Potatoes" 
        [set Qcanal item 7 wdM3-list  ]
        
        if current-crop = "Melon"
        [set Qcanal item 8 wdM3-list  ]
        
        if current-crop = "Carrots" 
        [set Qcanal item 9 wdM3-list  ]
        
        if current-crop = "Citrus" 
        [set Qcanal item 10 wdM3-list ]
        
        if current-crop = "fallow" 
        [set Qcanal item 11 wdM3-list ]
      
  
    set canal-water-used canal-water-used + Qcanal     
end

to UPDATE-PUMPING-RATES
 ; set well-water-used 0
 ; ask patches with [interior-node? = true and injection? = false and dry-well? = false]  ;;; MATCH COLOR TO ACTUAL CROP COLOR
 ; [
;if canal-water? = false[

;; rate of flow for groundwater is daily, so delta-t is used
  if current-crop = "Grains" 
    [set Qcanal item 0 wdM3-list ]   
  
  if current-crop = "Alfalfa"
    [set Qcanal item 1 wdM3-list  ]
  
  if current-crop = "Cotton" 
    [set Qcanal item 2 wdM3-list  ]
  
  if current-crop = "Onion" 
    [set Qcanal item 3 wdM3-list  ]
  
  if current-crop = "Truck" 
    [set Qcanal item 4 wdM3-list  ]
  
  
  if current-crop = "Deciduous" 
    [set Qcanal item 5 wdM3-list  ]
  
  if current-crop = "Grapes" 
    [set Qcanal item 6 wdM3-list ]
  
  if current-crop = "Potatoes" 
    [set Qcanal item 7 wdM3-list  ]
  
  if current-crop = "Melon"
    [set Qcanal item 8 wdM3-list  ]
  
  if current-crop = "Carrots" 
    [set Qcanal item 9 wdM3-list  ]
  
  if current-crop = "Citrus" 
    [set Qcanal item 10 wdM3-list  ]
  
  if current-crop = "fallow" 
    [set Qcanal item 11 wdM3-list ]
  
  ;; record total groundwater used in the period
  
  set well-water-used well-water-used + Qwell * delta-t
end

to UPDATE-CROP-COUNTS
  set m-list [current-crop] of patches with [farmer-here? = true]
  set n-list [ ]
  foreach m-list
  [
    if (? = "Citrus") [set n-list lput 0 n-list]
    if (? = "Carrots") [set n-list lput 1 n-list]
    if (? = "Melon") [set n-list lput 2 n-list]
    if (? = "Potatoes") [set n-list lput 3 n-list]
    if (? = "Grapes") [set n-list lput 4 n-list]
    if (? = "Deciduous") [set n-list lput 5 n-list]
    if (? = "Truck") [set n-list lput 6 n-list]
    if (? = "Onions") [set n-list lput 7 n-list]
    if (? = "Cotton") [set n-list lput 8 n-list]
    if (? = "Grains") [set n-list lput 9 n-list]
    if (? = "Alfalfa") [set n-list lput 10 n-list]

  ]
  set ready-to-histogram? true
end

to UPDATE-CROP-VIEW
  ask patches with [farmer-here? = true] ;[interior-node? = true and injection? = false]
  [
    
    if current-crop = "Citrus" [set pcolor 27]
    if current-crop = "Carrots" [set pcolor 25]
    if current-crop = "Melon" [set pcolor 15]
    if current-crop = "Potatoes" [set pcolor 37]
    if current-crop = "Grapes" [set pcolor 125]
    if current-crop = "Deciduous" [set pcolor 68]
    if current-crop = "Truck" [set pcolor 56]
    if current-crop = "Onions" [set pcolor 117]
    if current-crop = "Cotton" [set pcolor 39]
    if current-crop = "Grains" [set pcolor 32]
    if current-crop = "Alfalfa" [set pcolor 64]

  ] 
end


;to-report max-item [ #gmlist ]
;  let $max-value first #gmlist
;  let $max-item 0
;  let $item 1
;  foreach but-first #gmlist
;  [ if ? > $max-value [ set $max-value ? set $max-item $item]
;    set $item $item + 1
;  ]
;  report $max-item 
;end

to-report max-item [ #list ] 
   report position (max #list) #list 
end 

to SELL-CROPS
  ask patches with [farmer-here? and dry-well? = false]
  [
    set money money + item current-crop-index gmlist   
    if bankrupt? = true[
      copy  
    ]
    
  ]
  ;[set counter 0];; choose index number
end


to UPDATE-CROP-DISTRIBUTION
  set %-of-Grains (count patches with [current-crop = "Grains"] / count patches with [farmer-here? = true]) 
  set %-of-Alfalfa (count patches with [current-crop = "Alfalfa"] / count patches with [farmer-here? = true]) 
  set %-of-Cotton (count patches with [current-crop = "Cotton"] / count patches with [farmer-here? = true]) 
  set %-of-Onion (count patches with [current-crop = "Onion"] / count patches with [farmer-here? = true]) 
  set %-of-Truck (count patches with [current-crop = "Truck"] / count patches with [farmer-here? = true]) 
  set %-of-Deciduous (count patches with [current-crop = "Deciduous"] / count patches with [farmer-here? = true]) 
  set %-of-Grapes (count patches with [current-crop = "Grapes"] / count patches with [farmer-here? = true]) 
  set %-of-Potatoes (count patches with [current-crop = "Potatoes"] / count patches with [farmer-here? = true]) 
  set %-of-Melon (count patches with [current-crop = "Melon"] / count patches with [farmer-here? = true]) 
  set %-of-Carrots (count patches with [current-crop = "Carrots"] / count patches with [farmer-here? = true]) 
  set %-of-Citrus (count patches with [current-crop = "Citrus"] / count patches with [farmer-here? = true]) 
  set %-of-Fallow (count patches with [current-crop = "Fallow"] / count patches with [farmer-here? = true]) 
end


to UPDATE-CANAL-WATER-PRICE
  ; 2015 UCCE ORANGES
  ; "District water is delivered via canal to the farmat a cost of $114 per acare-foot or 9.50
  ;  per acre-inch. Water costs are highly variable among districts and in drought years water 
  ;  may increase to as high as $1,000 to $1,800 per acre-foot. This study assumes a year with
  ;  normal water costs.
  set current-canal-volume max-canal-flow * drought-index
  ;  set canal-water-price current-canal-volume * (-7.75 / 2370852) + 100;  set canal-water-price current-canal-volume * (-58 / 2370852) + 100
  
  ;  min price of 5 dollars per acre-inch
  ;  price is in acre inch, but flows set in m3
  
  ;; this equation automatically creates a line to sample from given values for max-canal-flow and max-canal-price
  if max-canal-flow > 0 [set canal-water-price current-canal-volume * -1 / (max-canal-flow / max-canal-price) + max-canal-price + 1]
end

to set-pumping-cost-ex-ante
  ;; Pumping costs are in meters cubed; are converted under CALCULATE-MARGIN
  set pumping-cost electricity-price * 9.81 * (Hground - H) + water-price;; pumping-cost units m^2/(s^2 * m^3)
  set H-ex-ante H
end

to set-pumping-costs-ex-post
  let mean-H (H + H-ex-ante) / 2
  set pumping-cost electricity-price * 9.81 * (Hground - mean-H) + water-price
end 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; DRILLER PROCEDURES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to UPDATE-MAX-WATER-DEPTH
  set max-depth-list lput max [Hground - H] of patches with [farmer-here? = true] max-depth-list
  set max-depth-value max max-depth-list 
end

to CHECK-IF-DRY-WELLS
ask patches with [farmer-here? = true]
[
ifelse (Hground - H) < well-depth [set dry-well? FALSE][set dry-well? TRUE set new-well-needed? TRUE]
]
end

to DRILL-A-WELL
  ask patches with [farmer-here? = TRUE and new-well-needed? = TRUE]
  [
    set well-depth max-depth-value + additional-drilling-depth                                 
    set new-well-needed? FALSE
    set total-wells-drilled total-wells-drilled + 1
    ;set money money - drilling-cost
    set mortgage? TRUE
    set payments-made 0
  ]
end

to REPAY-WELL
  if mortgage? = TRUE [set money money - drilling-cost / 7]
  set payments-made payments-made + 1
  if payments-made = 7 [set mortgage? false]
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; RECHARGE PROCEDURES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to UPDATE-INFLOW-TO-BASIN
  ask patches with [injection? = true]
  [
   set Qinjection inflow-to-basin * drought-index
  ] 
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; GENERATE NAMES AND WRITE FILES FOR BEHAVIORSPACE EXPERIMENTS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to output-results
set cost-differentials_name ""
set price-sensitivity_name ""
set mutation-rate_name ""
set starting-water-price_name ""
set electricity-price_name ""
set %-optimizers_name ""
set target-water-level_name ""
set canal-water-price_name ""
set drought-index_name ""
set prices_name ""
set inflow-to-basin_name ""


set price-sensitivity_name word price-sensitivity "_"
set mutation-rate_name word mutation-rate "_"
set starting-water-price_name word starting-water-price "_"
set electricity-price_name word electricity-price "_"
set %-optimizers_name word %-optimizers "_"
set target-water-level_name word target-water-level "_"
set canal-water-price_name word canal-water-price "_"
set drought-index_name word drought-index "_"
set prices_name word prices "_"
set inflow-to-basin_name word inflow-to-basin "_"

set filename (word cost-differentials_name price-sensitivity_name mutation-rate_name starting-water-price_name electricity-price_name %-optimizers_name target-water-level_name canal-water-price_name drought-index_name prices_name inflow-to-basin_name behaviorspace-run-number)

export-all-plots (word filename "_plots.txt")
export-interface (word filename "_interface.png")
export-view (word filename "_world.png")

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; READ CROP DATASET ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to read-crop-data
  file-close-all
  file-open "PARAMETERS.csv"
  let header file-read-line
  set fileList []
  
  if file-at-end? [ stop ]   ; stop reading if it reachs end of file  
  let i 0
  while [i < 7]
  [ set csv file-read-line
    set csv word csv ","  ; add semi-colon at the end of the line (csv) for loop termination
    
    let mylist []  ; list of values
    while [not empty? csv]
    [ let $x position "," csv
      let $item substring csv 0 $x  ; extract items one by one
      ;show $item
      carefully [set $item read-from-string $item][] ; if item is "number", convert to number
      set mylist lput $item mylist  ; append to mylist
      set csv substring csv ($x + 1) length csv  ; remove one item and comma at each step
    ]
    set i i + 1
    set fileList lput mylist fileList    ; append each list to fileList (fileList is list of lists)
  ]
  
  ; extract each sub-list 
  set cropnames-list item 0 fileList
  set acreage-list item 1 fileList
  set yield-list item 2 fileList
  set wdAF-list item 3 fileList
  set wdM3-list item 4 fileList
  set netprice-list item 5 fileList
  set tree-grower? item 6 fileList
  if tree-grower? = 0 [set tree-grower? false]
  if tree-grower? = 1 [set tree-grower? true]
  
  set acreage-list-pct map [? / (sum acreage-list)] acreage-list 
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; LORENZ AND GINI ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  

to update-lorenz-and-gini
  ask patches with [farmer-here? = true] [ 
    set wealth money
    if wealth < 0 [set wealth 0]
  ]
  let num-people count patches with [farmer-here? = true]
  let sorted-wealths sort [wealth] of patches with [farmer-here? = true]
  let total-wealth sum sorted-wealths
  let wealth-sum-so-far 0
  let index 0
  set gini-index-reserve 0
  set lorenz-points []
  repeat num-people [
    set wealth-sum-so-far (wealth-sum-so-far + item index sorted-wealths)
    if total-wealth = 0 [set total-wealth 0.00001]
    set lorenz-points lput ((wealth-sum-so-far / total-wealth) * 100) lorenz-points
    set index (index + 1)
    set gini-index-reserve
    gini-index-reserve +
    (index / num-people) -
    (wealth-sum-so-far / total-wealth)
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; GIS MAP STUFF ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to prepare-GIS-map
  ask patches with [gis:contained-by? self kern-shapefile = false]
[
set fixed-head? false
set no-flow? true
set interior-node? false
set pcolor black
]
end


to prepare-output
  if ticks > 100
  [
  set output-almonds []
  set output-almonds lput (count patches with [current-crop = "Almonds" and interior-node? = true and injection? = false] / count patches with [interior-node? = true and injection? = false]) output-almonds
  set output-cherries []
  set output-cherries lput (count patches with [current-crop = "Cherries" and interior-node? = true and injection? = false] / count patches with [interior-node? = true and injection? = false]) output-cherries
  set output-citrus []
  set output-citrus lput (count patches with [current-crop = "Citrus" and interior-node? = true and injection? = false] / count patches with [interior-node? = true and injection? = false]) output-citrus
  set output-grapes []
  set output-grapes lput (count patches with [current-crop = "Grapes" and interior-node? = true and injection? = false] / count patches with [interior-node? = true and injection? = false]) output-grapes
  set output-pistachios []
  set output-pistachios lput (count patches with [current-crop = "Pistachios" and interior-node? = true and injection? = false] / count patches with [interior-node? = true and injection? = false]) output-pistachios
  set output-tomatoes []
  set output-tomatoes lput (count patches with [current-crop = "Tomatoes" and interior-node? = true and injection? = false] / count patches with [interior-node? = true and injection? = false]) output-tomatoes
  set output-cotton []
  set output-cotton lput (count patches with [current-crop = "Cotton" and interior-node? = true and injection? = false] / count patches with [interior-node? = true and injection? = false]) output-cotton
  set output-alfalfa []
  set output-alfalfa lput (count patches with [current-crop = "Alfalfa" and interior-node? = true and injection? = false] / count patches with [interior-node? = true and injection? = false]) output-alfalfa
  set output-silageandforage []
  set output-silageandforage lput (count patches with [current-crop = "Silage and Forage" and interior-node? = true and injection? = false] / count patches with [interior-node? = true and injection? = false]) output-silageandforage
  set output-onions []
  set output-onions lput (count patches with [current-crop = "Onions" and interior-node? = true and injection? = false] / count patches with [interior-node? = true and injection? = false]) output-onions
  set output-bellpeppers []
  set output-bellpeppers lput (count patches with [current-crop = "Peppers Bell" and interior-node? = true and injection? = false] / count patches with [interior-node? = true and injection? = false]) output-bellpeppers
  set output-potatoes []
  set output-potatoes lput (count patches with [current-crop = "Potatoes" and interior-node? = true and injection? = false] / count patches with [interior-node? = true and injection? = false]) output-potatoes
  set output-fallow []
  set output-fallow lput (count patches with [current-crop = "fallow" and interior-node? = true and injection? = false] / count patches with [interior-node? = true and injection? = false]) output-fallow
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
580
52
1071
564
-1
-1
14.6
1
10
1
1
1
0
0
0
1
0
32
0
32
1
1
1
ticks
30.0

BUTTON
300
85
469
118
1-SETUP GW
setup-model
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
58
1520
215
1553
K
K
0.01
100
1
0.01
1
[m/day]
HORIZONTAL

SLIDER
826
1015
954
1048
delta-t
delta-t
1
365
365
1
1
[days]
HORIZONTAL

TEXTBOX
62
1502
286
1520
HYDRAULIC PARAMETERS
12
0.0
1

INPUTBOX
792
1115
866
1175
N
31
1
0
Number

INPUTBOX
882
1115
954
1175
M
31
1
0
Number

INPUTBOX
982
1155
1065
1220
delta
800
1
0
Number

TEXTBOX
606
1001
795
1020
SIMULATION PARAMETERS
12
0.0
1

CHOOSER
642
860
734
905
left-bc
left-bc
"no-flow" "fixed-head"
1

CHOOSER
918
860
1010
905
right-bc
right-bc
"no-flow" "fixed-head"
0

CHOOSER
826
860
918
905
top-bc
top-bc
"no-flow" "fixed-head"
0

CHOOSER
732
860
828
905
bottom-bc
bottom-bc
"no-flow" "fixed-head"
0

INPUTBOX
826
905
918
965
top-bc-head
100
1
0
Number

INPUTBOX
642
905
734
965
left-bc-head
450
1
0
Number

INPUTBOX
918
905
1010
965
right-bc-head
100
1
0
Number

INPUTBOX
732
905
829
965
bottom-bc-head
0
1
0
Number

MONITOR
792
1175
866
1220
Size X [km]
(delta * N) / 1000
17
1
11

MONITOR
882
1175
954
1220
Size Y [km]
(delta * M) / 1000
17
1
11

CHOOSER
850
1530
1116
1575
view
view
"crops" "wells" "K (values)" "K (contours)" "T (values)" "S (values)" "heads (values)" "heads (contours)" "areal recharge" "ET" "DRAIN conductance" "DRAIN flow [m3/d] and head values" "DRAIN flow [L/s] and head values" "RIV flow [m3/d] and head values" "RIV flow [L/s] and head values" "RIV flow [L/s] and drawdowns" "RIV storage and head values" "RIV bottom and head values" "RIV conductance" "patch numbering" "boundary conditions"
6

BUTTON
1116
1550
1225
1583
NIL
reset-view
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
712
1015
824
1060
solver
solver
"steady-state" "transient"
1

INPUTBOX
952
1015
1040
1075
initial-heads
10
1
0
Number

SLIDER
58
1552
215
1585
aquifer-thickness
aquifer-thickness
1
400
50
1
1
[m]
HORIZONTAL

CHOOSER
602
1015
714
1060
aquifer-type
aquifer-type
"confined" "unconfined"
0

BUTTON
302
254
471
287
5-RUN
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
216
1520
292
1585
S
1.0E-4
1
0
Number

BUTTON
282
1846
436
1879
place pumping well
place-pumping-well
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
282
1880
436
1913
place injection well
place-injection-well
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
674
27
1202
56
FlowLogo Platform — Arvin-Edison Water Storage District 
18
0.0
1

TEXTBOX
62
1828
425
1846
DISCHARGE AND RECHARGE WELLS
12
0.0
1

TEXTBOX
608
1852
790
1908
(input well discharge or injection rate, then press button and click on model window to add well or injection point)
11
0.0
1

TEXTBOX
818
1080
968
1098
MODEL DIMENSIONS
12
0.0
1

TEXTBOX
828
1100
919
1118
Number of cells
11
0.0
1

TEXTBOX
976
1506
1012
1524
VIEWS
12
0.0
1

INPUTBOX
292
1520
368
1585
Sy
0.3
1
0
Number

BUTTON
302
130
471
163
2-INITIAL HEADS
initialize-heads
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
872
1150
887
1168
X
11
0.0
1

TEXTBOX
986
1138
1094
1157
Size of cells [m]
11
0.0
1

BUTTON
434
1846
519
1879
hide patches
hide-stresses-patches
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
58
1846
131
1912
Qwell-input
1000
1
0
Number

INPUTBOX
132
1846
205
1912
Qinjection-input
500
1
0
Number

BUTTON
434
1880
519
1913
hide labels
hide-stresses-labels
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
518
1846
603
1912
clear all wells
clear-stresses
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
866
1001
918
1019
timestep
11
0.0
1

TEXTBOX
224
1586
285
1604
confined S
11
0.0
1

TEXTBOX
292
1584
366
1602
unconfined S
11
0.0
1

TEXTBOX
1348
1518
1456
1548
AREAL RECHARGE
12
0.0
1

INPUTBOX
1346
1534
1435
1600
areal-recharge
0.2
1
0
Number

TEXTBOX
1348
1604
1582
1622
Recharge units: [m3/m2*day] = [m/day]
11
0.0
1

TEXTBOX
62
1916
253
1944
Pumping/Injection units: [m3/day]\nScreen level: [m] above reference
11
0.0
1

INPUTBOX
202
1846
281
1912
screen-level-input
20
1
0
Number

TEXTBOX
60
1718
210
1748
BOUNDARY CONDITIONS
12
0.0
1

BUTTON
192
1736
407
1769
change this cell to: fixed head
change-to-fixed-head
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
58
1736
194
1802
fixed-head-value-new
400
1
0
Number

BUTTON
192
1770
407
1803
change this cell to: no-flow
change-to-no-flow
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1054
1648
1222
1681
export heads to file
write-output-heads
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
212
1636
346
1669
set K for this cell
set-K
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
58
1636
213
1696
K-input
10
1
0
Number

TEXTBOX
58
1618
344
1648
HYDRAULIC CONDUCTIVITY
12
0.0
1

BUTTON
1054
1744
1223
1779
NIL
reset-ticks
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
366
1520
619
1585
set these hydraulic parameters for all cells
setup-hydraulic-parameters
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
1800
1548
2000
1581
sine-recharge-multiplier?
sine-recharge-multiplier?
1
1
-1000

BUTTON
1634
1534
1784
1600
clear all recharge values
clear-recharge
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1434
1534
1634
1567
set recharge for all patches
set-recharge-all-patches
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1434
1568
1634
1601
set this value for a single patch
set-recharge-single-patch
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
798
1896
1104
1980
This multiplier is a sine wave of recharge centered around the month of June. Rain is zero during the summer, end of spring and the begining of autumn. This multiplier is applied to all the patches with recharge. Settings can be modified in the code
11
0.0
1

SWITCH
804
1856
988
1889
sine-well-multiplier?
sine-well-multiplier?
1
1
-1000

TEXTBOX
1350
1834
1601
1852
EVAPOTRANSPIRATION
12
0.0
1

SLIDER
1560
1850
1732
1883
ET-extinction-depth
ET-extinction-depth
0
10
5
1
1
NIL
HORIZONTAL

INPUTBOX
1350
1850
1480
1910
ET-land-surface-elevation
50
1
0
Number

BUTTON
1730
1850
1836
1883
add ET cell
set-ET-patch
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1836
1850
1940
1883
clear all ET cells
clear-ET-patches-all
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
1480
1850
1562
1910
ET-max-rate
0.0010
1
0
Number

TEXTBOX
1356
1910
1506
1928
ET units = [m/day]
11
0.0
1

BUTTON
1060
1837
1228
1870
import K grid file
import-SGS-grid\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
1350
1632
1525
1662
DRAIN/SPRING
12
0.0
1

INPUTBOX
1348
1648
1443
1708
DRAIN-elevation
40
1
0
Number

INPUTBOX
1444
1648
1557
1708
DRAIN-conductance
500
1
0
Number

TEXTBOX
1350
1708
1542
1726
Drain conductance units: [m2/day]
11
0.0
1

BUTTON
1556
1648
1670
1681
add drain cell
set-drain-patch
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
1506
1754
1604
1814
RIV-conductance
50
1
0
Number

INPUTBOX
1350
1754
1431
1814
RIV-elevation
50
1
0
Number

INPUTBOX
1430
1754
1509
1814
RIV-bottom
47
1
0
Number

TEXTBOX
1350
1736
1500
1754
LOSING/GAINING RIVER CONDITION
12
0.0
1

BUTTON
1604
1754
1706
1787
add RIV cell
set-RIV-patch
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1670
1648
1793
1681
clear all DRAIN cells
clear-DRAIN-patches-all
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1704
1754
1807
1787
clear all RIV cells
clear-RIV-patches-all
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1116
1518
1226
1551
NIL
refresh-view
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1054
1678
1222
1711
export world
export-world \"exampleGWABM\"
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1054
1710
1222
1743
import world
import-world \"exampleGWABM\"\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1062
1787
1221
1822
NIL
initialize-ABM
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
303
43
470
76
social-model
social-model
0
1
-1000

MONITOR
1179
58
1254
103
YEAR
YEAR
17
1
11

CHOOSER
688
1530
830
1575
view-farmers
view-farmers
"normal" "stress levels"
0

INPUTBOX
56
1984
141
2062
Qfixed-flux
1000
1
0
Number

BUTTON
142
2002
292
2037
place fixed flux cell
place-fixed-flux
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
58
1966
246
1989
FIXED-FLUX
12
0.0
1

SLIDER
2347
909
2520
942
number-of-drillers
number-of-drillers
0
100
0
1
1
NIL
HORIZONTAL

BUTTON
303
210
473
245
4-SETUP SOCIAL
SETUP-SOCIAL
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
19
40
223
73
%-optimizers
%-optimizers
0
1
0
0.01
1
NIL
HORIZONTAL

PLOT
1378
712
2083
862
Reservoir Level
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" "set-plot-y-range 0 500"
PENS
"H" 1.0 0 -16777216 true "" "if model-ready? = true [plot mean [H] of patches with [interior-node? = true]]"
"Hinitial" 1.0 0 -7500403 true "" "plot 500"
"Target Water Level" 1.0 0 -2674135 true "" "if model-ready? = true [plot target-water-level * mean [Hground] of patches with [interior-node? = true]]"
"mean well depth" 1.0 2 -14070903 true "" "if model-ready? = true [plot mean [Hground - well-depth] of patches with [interior-node? = true]]"

PLOT
1838
552
2077
704
Water Prices
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Ground Water" 1.0 0 -16777216 true "" "plot water-price"
"Canal Water" 1.0 0 -2674135 true "" "plot canal-water-price"

SLIDER
27
798
192
831
target-water-level
target-water-level
0
1
0.8
.01
1
NIL
HORIZONTAL

SLIDER
19
79
223
112
remaining-%-copycats
remaining-%-copycats
0
1
1
.01
1
NIL
HORIZONTAL

PLOT
1378
553
1618
703
Farmer Wealth Distribution
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" "set-histogram-num-bars 20\nset-plot-x-range 0 (max [money] of patches + .0001)\nset-plot-pen-interval ((max [money] of patches) / 50 + .0001)"
PENS
"default" 1.0 1 -16777216 true "" "histogram [round money] of patches"

SWITCH
27
707
192
740
prices
prices
1
1
-1000

SLIDER
28
840
190
873
price-sensitivity
price-sensitivity
0
500
10
10
1
NIL
HORIZONTAL

SLIDER
269
875
473
908
drilling-cost
drilling-cost
0
1000000
1000000
10000
1
[US$]
HORIZONTAL

SLIDER
19
122
222
155
mutation-rate
mutation-rate
0
1
0.05
.01
1
NIL
HORIZONTAL

MONITOR
1748
294
1825
339
NIL
water-price
0
1
11

PLOT
1370
395
1617
545
Crop Distribution
NIL
NIL
0.0
10.0
0.0
400.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" "if ready-to-histogram? = true [histogram n-list]"

PLOT
1623
395
1831
545
Lorenz Curve
Pop %
Wealth %
0.0
100.0
0.0
100.0
false
false
"" ""
PENS
"lorenz" 1.0 0 -16777216 true "" "if model-ready? = true\n[\nplot-pen-reset\nset-plot-pen-interval 100 / count patches with [interior-node? = true and injection? = false]\nplot 0\nforeach lorenz-points plot\n]"
"equal" 100.0 0 -7500403 true "plot 0\nplot 100" ""

PLOT
1368
44
2092
407
% of crops
NIL
NIL
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"Grains" 1.0 0 -12440034 true "" "if model-ready? = true [plot count patches with [current-crop = \"Grains\" and farmer-here? = true] / count patches with [farmer-here? = true]]"
"Alfalfa" 1.0 0 -14439633 true "" "if model-ready? = true [plot count patches with [current-crop = \"Alfalfa\" and farmer-here? = true] / count patches with [farmer-here? = true]]"
"Cotton" 1.0 0 -1318182 true "" "if model-ready? = true [plot count patches with [current-crop = \"Cotton\" and farmer-here? = true] / count patches with [farmer-here? = true]]"
"Onion" 1.0 0 -5204280 true "" "if model-ready? = true [plot count patches with [current-crop = \"Onion\" and farmer-here? = true] / count patches with [farmer-here? = true]]"
"Truck" 1.0 0 -8732573 true "" "if model-ready? = true [plot count patches with [current-crop = \"Truck\" and farmer-here? = true] / count patches with [farmer-here? = true]]"
"Deciduous" 1.0 0 -5509967 true "" "if model-ready? = true [plot count patches with [current-crop = \"Deciduous\" and farmer-here? = true] / count patches with [farmer-here? = true]]"
"Grapes" 1.0 0 -5825686 true "" "if model-ready? = true [plot count patches with [current-crop = \"Grapes\" and farmer-here? = true] / count patches with [farmer-here? = true]]"
"Potatoes" 1.0 0 -3889007 true "" "if model-ready? = true [plot count patches with [current-crop = \"Potatoes\" and farmer-here? = true] / count patches with [farmer-here? = true]]"
"Melon" 1.0 0 -2674135 true "" "if model-ready? = true [plot count patches with [current-crop = \"Melon\" and farmer-here? = true] / count patches with [farmer-here? = true]]"
"Carrots" 1.0 0 -955883 true "" "if model-ready? = true [plot count patches with [current-crop = \"Carrots\" and farmer-here? = true] / count patches with [farmer-here? = true]]"
"Citrus" 1.0 0 -612749 true "" "if model-ready? = true [plot count patches with [current-crop = \"Citrus\" and farmer-here? = true] / count patches with [farmer-here? = true]]"

BUTTON
2347
1093
2517
1127
2-LOAD KERN OVERLAY
set kern-shapefile gis:load-dataset \"kern2.shp\"\nset gis-envelope gis:envelope-of kern-shapefile\ngis:set-world-envelope-ds gis-envelope\ngis:set-drawing-color black\ngis:draw kern-shapefile 2\nask patches with [gis:contained-by? self kern-shapefile = false]\n[\nset fixed-head? false\nset no-flow? true\nset pcolor black\n]\nimport-drawing \"kern mask.png\"
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
1624
553
1831
701
Total Farmer Wealth
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot sum [wealth] of patches with [interior-node? = true and injection? = false]"

SWITCH
272
702
372
735
drilling?
drilling?
0
1
-1000

MONITOR
225
410
326
455
mean well depth
mean [well-depth] of patches with [farmer-here? = true]
1
1
11

PLOT
1838
395
2078
547
wells drilled
NIL
number of wells
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot total-wells-drilled"

SLIDER
20
162
222
195
drought-index
drought-index
0
1
0.75
0.05
1
NIL
HORIZONTAL

SLIDER
20
202
220
235
electricity-price
electricity-price
0
1
0.15
0.01
1
[$/kWh]
HORIZONTAL

MONITOR
20
247
220
292
canal water price [$/ac-ft]
canal-water-price
0
1
11

TEXTBOX
74
16
191
35
PARAMETERS
18
0.0
1

TEXTBOX
337
17
469
62
SIMULATION
18
0.0
1

TEXTBOX
29
645
219
719
DYNAMIC SCARCITY PRICING
18
0.0
1

BUTTON
302
170
474
205
3-SETUP SPREADING BASINS
ask patch 15 15\n[ \nset gw-bank-here? true\nset farmer-here? false\nset pcolor red \n]\n\nask patches with [interior-node? = true]\n[\nset H 400\n]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
302
299
474
332
water-level-labels?
water-level-labels?
0
1
-1000

SLIDER
20
365
222
398
inflow-to-basin
inflow-to-basin
0
10000
4500
50
1
[m3/day]
HORIZONTAL

MONITOR
124
408
219
453
Total outflows
sum [Qwell] of patches with [interior-node? = true and injection? = false]
0
1
11

MONITOR
25
407
118
452
Total inflows
sum [Qinjection] of patches with [injection? = true]
0
1
11

SLIDER
20
329
219
362
max-canal-price
max-canal-price
0
1800
1800
1
1
NIL
HORIZONTAL

SLIDER
27
753
196
786
starting-water-price
starting-water-price
0
800
5
5
1
NIL
HORIZONTAL

SLIDER
29
885
191
918
median-starting-farmer-money
median-starting-farmer-money
0
1000000
0
10000
1
NIL
HORIZONTAL

SWITCH
24
524
316
557
starting-wealth-distribution-uniform?
starting-wealth-distribution-uniform?
1
1
-1000

SLIDER
30
928
191
961
standard-deviation-size
standard-deviation-size
0
100
50
1
1
NIL
HORIZONTAL

PLOT
1378
869
2083
1019
Average Pumping Cost and Price Plus Pumping Cost per Acre-foot
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Mean Pumping Costs" 1.0 0 -16777216 true "" "if plots-ready? = true\n[plot mean [pumping-cost] of patches with [interior-node? = true]]"
"Max Pumping Costs" 1.0 0 -7500403 true "" "if plots-ready? = true\n[plot max [pumping-cost ] of patches with [interior-node? = true]]"
"Min Pumping Costs" 1.0 0 -2674135 true "" "if plots-ready? = true\n[plot min [pumping-cost ] of patches with [interior-node? = true]]"
"pen-3" 1.0 0 -955883 true "" ""

SLIDER
21
293
219
326
max-canal-flow
max-canal-flow
0
4000000
3088000
1000
1
NIL
HORIZONTAL

PLOT
1378
1027
2085
1177
Water Use
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Canal" 1.0 0 -16777216 true "" "plot canal-water-used"
"Well" 1.0 0 -7500403 true "" "plot well-water-used"
"Total" 1.0 0 -2674135 true "" "plot well-water-used + canal-water-used"

PLOT
1375
1208
2127
1358
Dry Wells
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ";if model-ready? = true \n;[plot count patches with [interior-node? and dry-well?]]"

PLOT
2097
399
2454
621
Bankruptcies
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count-bankrupt"

PLOT
2097
619
2455
832
Bankruptcies per Tick
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 2 -16777216 true "" "plot count patches with [bankrupt? = true]"

SLIDER
272
748
475
781
initial-well-depth
initial-well-depth
0
500
200
10
1
[m]
HORIZONTAL

SLIDER
268
834
475
867
waiting-time-for-well
waiting-time-for-well
0
3
1
1
1
[years]
HORIZONTAL

PLOT
2099
54
2519
399
Strategies
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Optimization" 1.0 0 -16777216 true "" "plot count patches with [interior-node? = true and optimizer? = true]"
"Copy Cat" 1.0 0 -7500403 true "" "plot count patches with [interior-node? = true and copycat? = true]"
"Adaptive Only" 1.0 0 -2674135 true "" "plot count patches with [interior-node? = true and optimizer? = false and copycat? = false]"

SWITCH
24
564
264
597
change-strategy-if-bankrupt?
change-strategy-if-bankrupt?
0
1
-1000

MONITOR
137
467
254
512
farm size [acres]
(delta * delta) / 10000 * 2.47\n;;average operation size is 158 acres considering ~112,617 irrigated acres
1
1
11

MONITOR
22
465
129
510
number of fams
count patches with [farmer-here? = true]\n;;AEWSD water management plan says ~711 farms/parcels
0
1
11

MONITOR
262
467
471
512
district agricultural area [acres]
count patches with [farmer-here? = true] * (delta * delta) / 10000 * 2.47\n;; AEWSD services ca. 116000 acres
0
1
11

MONITOR
464
384
604
429
KD WD water level [m]
left-bc-head
0
1
11

MONITOR
859
763
1018
808
WRM WSD water level [m]
bottom-bc-head
0
1
11

TEXTBOX
273
642
435
692
DRILLING DYNAMICS
16
0.0
1

INPUTBOX
1052
858
1207
918
buffer-zone
3
1
0
Number

SLIDER
269
790
502
823
additional-drilling-depth
additional-drilling-depth
0
100
50
1
1
[m]
HORIZONTAL

@#$#@#$#@
## THINGS TO DO

Acre to Hectare
Consolidate gm variables
Add slider for cost differentials


## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

Farmer wealth distriution is either uniform ( starting-wealth-distribution-uniform = true) or standard-normal ( starting-wealth-distribution-uniform = False ). The Uniform
distribution includes the range (0, 2 * median-starting-farmer-wealth). 

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

house bungalow
false
0
Rectangle -7500403 true true 210 75 225 255
Rectangle -7500403 true true 90 135 210 255
Rectangle -16777216 true false 165 195 195 255
Line -16777216 false 210 135 210 255
Rectangle -16777216 true false 105 202 135 240
Polygon -7500403 true true 225 150 75 150 150 75
Line -16777216 false 75 150 225 150
Line -16777216 false 195 120 225 150
Polygon -16777216 false false 165 195 150 195 180 165 210 195
Rectangle -16777216 true false 135 105 165 135

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

person construction
false
0
Rectangle -7500403 true true 123 76 176 95
Polygon -1 true false 105 90 60 195 90 210 115 162 184 163 210 210 240 195 195 90
Polygon -13345367 true false 180 195 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285
Circle -7500403 true true 110 5 80
Line -16777216 false 148 143 150 196
Rectangle -16777216 true false 116 186 182 198
Circle -1 true false 152 143 9
Circle -1 true false 152 166 9
Rectangle -16777216 true false 179 164 183 186
Polygon -955883 true false 180 90 195 90 195 165 195 195 150 195 150 120 180 90
Polygon -955883 true false 120 90 105 90 105 165 105 195 150 195 150 120 120 90
Rectangle -16777216 true false 135 114 150 120
Rectangle -16777216 true false 135 144 150 150
Rectangle -16777216 true false 135 174 150 180
Polygon -955883 true false 105 42 111 16 128 2 149 0 178 6 190 18 192 28 220 29 216 34 201 39 167 35
Polygon -6459832 true false 54 253 54 238 219 73 227 78
Polygon -16777216 true false 15 285 15 255 30 225 45 225 75 255 75 270 45 285

person farmer
false
0
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Polygon -1 true false 60 195 90 210 114 154 120 195 180 195 187 157 210 210 240 195 195 90 165 90 150 105 150 150 135 90 105 90
Circle -7500403 true true 110 5 80
Rectangle -7500403 true true 127 79 172 94
Polygon -13345367 true false 120 90 120 180 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 180 90 172 89 165 135 135 135 127 90
Polygon -6459832 true false 116 4 113 21 71 33 71 40 109 48 117 34 144 27 180 26 188 36 224 23 222 14 178 16 167 0
Line -16777216 false 225 90 270 90
Line -16777216 false 225 15 225 90
Line -16777216 false 270 15 270 90
Line -16777216 false 247 15 247 90
Rectangle -6459832 true false 240 90 255 300

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="false">
    <setup>setup-model
set kern-shapefile gis:load-dataset "kern2.shp"
set gis-envelope gis:envelope-of kern-shapefile
gis:set-world-envelope-ds gis-envelope
gis:set-drawing-color black
gis:draw kern-shapefile 2
ask patches with [gis:contained-by? self kern-shapefile = false]
[
set fixed-head? false
set no-flow? true
set pcolor black
]
import-drawing "kern mask.png"
initialize-heads
ask patch 35 20 
[ 
set pcolor magenta 
set fixed-flux? true
set injection? true
set Qinjection inflow-to-basin
]

ask patch 21 15
[ 
set pcolor magenta 
set fixed-flux? true
set injection? true
set Qinjection inflow-to-basin
]

ask patches with [interior-node? = true]
[
set H 400
]
SETUP-SOCIAL</setup>
    <go>go</go>
    <final>output-results</final>
    <timeLimit steps="1000"/>
    <metric>mean output-almonds</metric>
    <metric>mean output-almonds</metric>
    <metric>mean output-cherries</metric>
    <metric>mean output-citrus</metric>
    <metric>mean output-grapes</metric>
    <metric>mean output-pistachios</metric>
    <metric>mean output-tomatoes</metric>
    <metric>mean output-cotton</metric>
    <metric>mean output-alfalfa</metric>
    <metric>mean output-silageandforage</metric>
    <metric>mean output-onions</metric>
    <metric>mean output-bellpeppers</metric>
    <metric>mean output-potatoes</metric>
    <metric>mean output-fallow</metric>
    <enumeratedValueSet variable="cost-differentials">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="price-sensitivity">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="starting-water-price">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="electricity-price">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%-optimizers">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="target-water-level">
      <value value="0.75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="canal-water-price">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drought-index">
      <value value="0.2"/>
      <value value="0.25"/>
      <value value="0.3"/>
      <value value="0.35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prices">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inflow-to-basin">
      <value value="1250"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="top-bc">
      <value value="&quot;no-flow&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RIV-bottom">
      <value value="47"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="social-model">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cost-differentials">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aquifer-thickness">
      <value value="400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RIV-conductance">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="price-sensitivity">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DRAIN-conductance">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ET-max-rate">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-rate">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bottom-bc">
      <value value="&quot;no-flow&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="starting-water-price">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="right-bc-head">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fixed-head-value-new">
      <value value="400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="well-cost">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="electricity-price">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Qfixed-flux">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drilling-cost">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="solver">
      <value value="&quot;transient&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sine-recharge-multiplier?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-heads">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-t">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="right-bc">
      <value value="&quot;no-flow&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="left-bc-head">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ET-land-surface-elevation">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DRAIN-elevation">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="M">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aquifer-type">
      <value value="&quot;confined&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%-optimizers">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sine-well-multiplier?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="target-water-level">
      <value value="0.74"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="left-bc">
      <value value="&quot;no-flow&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="view-farmers">
      <value value="&quot;normal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="canal-water-price">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K-input">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-level-labels?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Sy">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="remaining-%-copycats">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="screen-level-input">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drought-index">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Qwell-input">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prices">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Qinjection-input">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="areal-recharge">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ET-extinction-depth">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inflow-to-basin">
      <value value="1250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="view">
      <value value="&quot;heads (values)&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RIV-elevation">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-drillers">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bottom-bc-head">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drilling?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="top-bc-head">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="top-bc">
      <value value="&quot;no-flow&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RIV-bottom">
      <value value="47"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="social-model">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cost-differentials">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aquifer-thickness">
      <value value="400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RIV-conductance">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="price-sensitivity">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DRAIN-conductance">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ET-max-rate">
      <value value="0.0010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-rate">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bottom-bc">
      <value value="&quot;no-flow&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="starting-water-price">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="right-bc-head">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fixed-head-value-new">
      <value value="400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="well-cost">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="electricity-price">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Qfixed-flux">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drilling-cost">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="solver">
      <value value="&quot;transient&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sine-recharge-multiplier?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-heads">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-t">
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="right-bc">
      <value value="&quot;no-flow&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="left-bc-head">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ET-land-surface-elevation">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DRAIN-elevation">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="M">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aquifer-type">
      <value value="&quot;confined&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%-optimizers">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sine-well-multiplier?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="target-water-level">
      <value value="0.75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="left-bc">
      <value value="&quot;no-flow&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="view-farmers">
      <value value="&quot;normal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="canal-water-price">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K-input">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-level-labels?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Sy">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="remaining-%-copycats">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="screen-level-input">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drought-index">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Qwell-input">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prices">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Qinjection-input">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="areal-recharge">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ET-extinction-depth">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inflow-to-basin">
      <value value="1250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="view">
      <value value="&quot;heads (values)&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RIV-elevation">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-drillers">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bottom-bc-head">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drilling?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="top-bc-head">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
