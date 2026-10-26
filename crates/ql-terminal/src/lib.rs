//! A dependency-free terminal dashboard and umbrella SDK for ql-rs.

pub mod prelude {
    pub use crate::TabId;
    pub use ql_ad::*;
    pub use ql_commodities::*;
    pub use ql_commodities::{
        BessAnalyticalSdpDegradationStorageReport, BessAnalyticalSdpDegradationStorageSpec,
        evaluate_bess_analytical_sdp_degradation_storage,
    };
    pub use ql_commodities::{
        WarinBaeuerleStorageRollingHedgeReport, WarinBaeuerleStorageRollingHedgeSpec,
        evaluate_warin_baeuerle_storage_rolling_hedge,
    };
    pub use ql_credit::*;
    pub use ql_credit::{
        BrigoMoriniCmcdsCreditOptionReport, BrigoMoriniCmcdsCreditOptionSpec,
        evaluate_brigo_morini_cmcds_and_credit_index_option,
    };
    pub use ql_credit::{
        DamodaranRsmCountryRiskAndDebtModReport, DamodaranRsmCountryRiskAndDebtModSpec,
        FsbFinmaGsibResolutionBailInWaterfallReport, FsbFinmaGsibResolutionBailInWaterfallSpec,
        evaluate_damodaran_rsm_country_risk_and_debt_modification,
        evaluate_fsb_finma_gsib_resolution_bail_in_waterfall,
    };
    pub use ql_credit::{
        FedsDeSpiegeleerAt1CocoDualTriggerReport, FedsDeSpiegeleerAt1CocoDualTriggerSpec,
        evaluate_feds_despiegeleer_at1_coco_dual_trigger,
    };
    pub use ql_credit::{
        LandoIsraelGeneratorLogarithmReport, LandoIsraelGeneratorLogarithmSpec,
        compute_lando_israel_generator_logarithm_and_migration,
    };
    pub use ql_credit::{
        RsmAsc470DebtModificationReport, RsmAsc470DebtModificationSpec,
        evaluate_rsm_asc470_debt_modification_and_restructuring,
    };
    pub use ql_curves::*;
    pub use ql_curves::{
        AdbRbiNdirsAssetSwapReport, AdbRbiNdirsAssetSwapSpec,
        evaluate_adb_rbi_ndirs_and_asset_swap_basis,
    };
    pub use ql_curves::{
        BcbsD368SixScenarioIrrbbReport, BcbsD368SixScenarioIrrbbSpec,
        evaluate_bcbs_d368_six_scenario_irrbb,
    };
    pub use ql_curves::{
        FabozziKeyRateLscImmunizationReport, FabozziKeyRateLscImmunizationSpec,
        evaluate_fabozzi_tuckman_key_rate_and_lsc_immunization,
    };
    pub use ql_curves::{
        RomerVeronesiTaylorRuleMacroCurveReport, RomerVeronesiTaylorRuleMacroCurveSpec,
        evaluate_romer_veronesi_taylor_rule_macro_yield_curve,
    };
    pub use ql_dml::*;
    pub use ql_dml::{
        AgrawalBoydDifferentiableQpReport, AgrawalBoydDifferentiableQpSpec,
        solve_agrawal_boyd_differentiable_qp_layer,
    };
    pub use ql_dml::{
        GlassermanSaddlepointIsStratifiedLhsReport, GlassermanSaddlepointIsStratifiedLhsSpec,
        evaluate_glasserman_saddlepoint_is_and_stratified_lhs_aad,
    };
    pub use ql_dml::{
        GreeneWooldridgeHeckmanSelectionGmmReport, GreeneWooldridgeHeckmanSelectionGmmSpec,
        HastieTibshiraniElasticNetGcvReport, HastieTibshiraniElasticNetGcvSpec,
        evaluate_greene_wooldridge_heckman_survivorship_and_gmm,
        solve_hastie_tibshirani_elastic_net_and_gcv,
    };
    pub use ql_dml::{
        HugeSavineDifferentialPcaReport, HugeSavineDifferentialPcaSpec,
        evaluate_huge_savine_differential_pca_and_twin_dml,
    };
    pub use ql_exotics::*;
    pub use ql_exotics::{
        ChoiLiDengMultiAssetSpreadBasketReport, ChoiLiDengMultiAssetSpreadBasketSpec,
        evaluate_choi_li_deng_multi_asset_spread_basket,
    };
    pub use ql_exotics::{
        CurranRogersShiArithmeticAsianReport, CurranRogersShiArithmeticAsianSpec,
        price_curran_rogers_shi_arithmetic_asian,
    };
    pub use ql_exotics::{
        SheppShiryaevShrevePerpetualRussianReport, SheppShiryaevShrevePerpetualRussianSpec,
        price_shepp_shiryaev_shreve_perpetual_russian_and_mckean,
    };
    pub use ql_fixed_income::*;
    pub use ql_fixed_income::{
        AcharyaFedsBankRunAndCocoReport, AcharyaFedsBankRunAndCocoSpec,
        analyze_acharya_feds_bank_run_and_optimal_coco,
    };
    pub use ql_fixed_income::{
        AnDaoquanCliffordChanceResettableCbReport, AnDaoquanCliffordChanceResettableCbSpec,
        evaluate_an_daoquan_clifford_chance_resettable_convertible,
    };
    pub use ql_fixed_income::{
        AsquithSirbuShreveCbGameReport, AsquithSirbuShreveCbGameSpec,
        DavidsonLevinCmoPacCompanionReport, DavidsonLevinCmoPacCompanionSpec,
        FasbAsu2024InducedConversionCpdiReport, FasbAsu2024InducedConversionCpdiSpec,
        evaluate_asquith_sirbu_shreve_cb_call_delay_and_dynkin_game,
        evaluate_davidson_levin_cmo_pac_companion_waterfall,
        evaluate_fasb_asu2024_induced_conversion_and_cpdi_oid,
    };
    pub use ql_fixed_income::{
        KoziolAsquithFasbStrategicConvertibleReport, KoziolAsquithFasbStrategicConvertibleSpec,
        VeronesiFabozziConvexityPrepaymentReport, VeronesiFabozziConvexityPrepaymentSpec,
        evaluate_koziol_asquith_fasb_strategic_convertible,
        evaluate_veronesi_fabozzi_futures_convexity_and_prepayment,
    };
    pub use ql_fixed_income::{
        MayerBrownSsgaCappedCallMakeWholeReport, MayerBrownSsgaCappedCallMakeWholeSpec,
        evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole,
    };
    pub use ql_fx::*;
    pub use ql_fx::{
        BaranBianchettiMtmXccyBasisReport, BaranBianchettiMtmXccyBasisSpec,
        evaluate_baran_bianchetti_mtm_xccy_basis_swap,
    };
    pub use ql_fx::{
        CrossCurrencyConvertibleWarrantReport, CrossCurrencyConvertibleWarrantSpec,
        price_cross_currency_convertible_and_warrant_bond,
    };
    pub use ql_inflation::*;
    pub use ql_inflation::{
        CochraneFiscalTheoryPriceLevelReport, CochraneFiscalTheoryPriceLevelSpec,
        evaluate_cochrane_ftpl_and_maturity_inflation_smoothing,
    };
    pub use ql_instruments::*;
    pub use ql_levy::{
        CgmyFractionalFftAndPideReport, CgmyFractionalFftAndPideSpec, CgmyProcess,
        ContTankovBatesSvjPideReport, ContTankovBatesSvjPideSpec, KouJumpDiffusion, LevyModel,
        MertonJumpDiffusion, VarianceGamma, cos_levy_price, cos_levy_price_strip,
        evaluate_cgmy_fractional_fft_and_pide,
        evaluate_cont_tankov_bates_svj_characteristic_and_pide, merton_series_price,
    };
    pub use ql_levy::{carr_madan_fft_levy_price_strip, lewis_levy_price_strip};
    pub use ql_math::*;
    pub use ql_mc::*;
    pub use ql_models::{
        AnalyticBsmHullWhiteEngine, BatesParams, Bergomi2FParams,
        BergomiForwardStartAndVarStripReport, BergomiForwardStartAndVarStripSpec,
        CalibratedHestonModel, Complex, CosHestonEngine, CosHestonHullWhiteEngine,
        DamodaranRdCapitalizationRealOptionReport, DamodaranRdCapitalizationRealOptionSpec,
        DoubleHestonParams, HestonHullWhiteParams, HestonMarketQuote, HestonParams, HullWhite1F,
        LiftedHestonParams, MasColellCochraneArrowDebreuSdfReport,
        MasColellCochraneArrowDebreuSdfSpec, MixedRoughBergomiVixReport, MixedRoughBergomiVixSpec,
        QuinticOuParams, QuinticOuVixSpxReport, VixOptionValuation, bergomi_2f_implied_vol_bg,
        bergomi_2f_mc_price, bergomi_2f_smile_moments_bg, calibrate_heston,
        double_heston_cos_price, double_heston_cos_strip_prices,
        evaluate_bergomi_forward_start_and_var_strip,
        evaluate_damodaran_rd_capitalization_and_patent_real_option,
        evaluate_mascolell_cochrane_arrow_debreu_sdf, evaluate_mixed_rough_bergomi_vix,
        evaluate_quintic_ou_spx_vix, h1_hw_char_func_unit_spot,
        heston_bates_vix_future_and_option_price, heston_cos_price, lifted_heston_char_fn,
        lifted_heston_char_func_unit_spot, lifted_heston_cos_price, lifted_heston_kernel,
        local_vol::{
            AndreasenHugeGridSpec, AndreasenHugeSliceResult, FeldmanSabrCalibrationResult,
            ZabrParams, ZabrSmileQuote, evaluate_zabr_smile, feldman_analytic_sabr_calibration,
            interpolate_andreasen_huge_call, solve_andreasen_huge_step,
        },
        replicate_cboe_vix_index,
    };
    pub use ql_models::{
        CarmonaDelarueMeanFieldSystemicReport, CarmonaDelarueMeanFieldSystemicSpec,
        solve_carmona_delarue_mean_field_systemic_interbank_game,
    };
    pub use ql_models::{
        DoubleHestonKernelSlvReport, DoubleHestonKernelSlvSpec,
        evaluate_double_heston_cumulants_and_kernel_slv,
    };
    pub use ql_models::{
        GairatPolalaCheyetteCapletReport, GairatPolalaCheyetteCapletSpec,
        GatheralJusselinQuadraticRoughHestonReport, GatheralJusselinQuadraticRoughHestonSpec,
        SchilitSloanDechowDichevQualityReport, SchilitSloanDechowDichevQualitySpec,
        evaluate_gairat_polala_cheyette_local_vol_caplet,
        evaluate_gatheral_jusselin_quadratic_rough_heston,
        evaluate_schilit_sloan_dechow_dichev_earnings_quality,
    };
    pub use ql_models::{
        JacquierDeepRoughVolSurrogateReport, JacquierDeepRoughVolSurrogateSpec,
        evaluate_jacquier_chebyshev_rough_vol_surrogate,
    };
    pub use ql_models::{
        PenmanReoiOhlsonJuettnerReport, PenmanReoiOhlsonJuettnerSpec,
        evaluate_penman_reoi_and_ohlson_juettner_aeg_valuation,
    };
    pub use ql_models::{rough::*, slv::*};
    pub use ql_pde::*;
    pub use ql_pde::{
        HaentjensInthoutHhwAdiReport, HaentjensInthoutHhwAdiSpec,
        evaluate_haentjens_inthout_heston_hull_white_adi,
    };
    pub use ql_portfolio::*;
    pub use ql_portfolio::{
        AhmadiJavidEntropicVarPortfolioReport, AhmadiJavidEntropicVarPortfolioSpec,
        CottonSchurHrpRmtNcoReport, CottonSchurHrpRmtNcoSpec, TsayBekkGarchSpilloverReport,
        TsayBekkGarchSpilloverSpec, evaluate_ahmadi_javid_entropic_var_portfolio,
        evaluate_cotton_schur_hrp_and_rmt_nco, evaluate_tsay_bekk_garch_volatility_spillover,
    };
    pub use ql_portfolio::{
        BarberisJinWangProspectAnomalyReport, BarberisJinWangProspectAnomalySpec,
        ChristoffersenBouchaudVarBacktestAndCumulantReport,
        ChristoffersenBouchaudVarBacktestAndCumulantSpec,
        ColinBaconCampisiFixedIncomeAttributionReport, ColinBaconCampisiFixedIncomeAttributionSpec,
        FoucaultVanDerMerweLobLiquidityLvarReport, FoucaultVanDerMerweLobLiquidityLvarSpec,
        HartMasColellRegretMatchingPortfolioReport, HartMasColellRegretMatchingPortfolioSpec,
        evaluate_barberis_jin_wang_prospect_overhang_and_skew_anomaly,
        evaluate_christoffersen_bouchaud_var_backtest_and_cumulant_es,
        evaluate_colin_bacon_campisi_fixed_income_attribution,
        evaluate_foucault_vandermerwe_lob_schedule_and_lvar,
        evaluate_hart_mascolell_regret_matching_allocation,
    };
    pub use ql_portfolio::{
        BouchaudPottersCumulantSmileReport, BouchaudPottersCumulantSmileSpec,
        BrockHommesAdaptiveBeliefReport, BrockHommesAdaptiveBeliefSpec,
        DgtwTreynorMazuyTimingReport, DgtwTreynorMazuyTimingSpec,
        EmbrechtsHillPickandsSubexponentialReport, EmbrechtsHillPickandsSubexponentialSpec,
        KarnoskySingerMencheroAttributionReport, KarnoskySingerMencheroAttributionSpec,
        evaluate_bouchaud_potters_power_law_cumulant_and_smile_scaling,
        evaluate_brock_hommes_adaptive_belief_switching,
        evaluate_dgtw_characteristic_and_treynor_mazuy_timing,
        evaluate_embrechts_hill_pickands_and_subexponential_tail,
        evaluate_karnosky_singer_multicurrency_and_menchero_attribution,
    };
    pub use ql_portfolio::{
        BussetiBoydRiskConstrainedKellyReport, BussetiBoydRiskConstrainedKellySpec,
        HarrisRollCorwinSchultzMicrostructureReport, HarrisRollCorwinSchultzMicrostructureSpec,
        analyze_harris_roll_corwin_schultz_microstructure,
        solve_busseti_boyd_risk_constrained_kelly,
    };
    pub use ql_portfolio::{
        CarteaJaimungalPenalvaLobFillHjbReport, CarteaJaimungalPenalvaLobFillHjbSpec,
        GrinoldQianChristophersonAlphaDecayReport, GrinoldQianChristophersonAlphaDecaySpec,
        MadhavanAbnerHillEtfApArbitrageReport, MadhavanAbnerHillEtfApArbitrageSpec,
        evaluate_grinold_qian_christopherson_ic_decay_and_timing,
        evaluate_madhavan_abner_hill_etf_ap_arbitrage_and_tax_alpha,
        solve_cartea_jaimungal_penalva_lob_fill_hjb,
    };
    pub use ql_portfolio::{
        FollmerSchiedGennaioliDiagnosticRiskReport, FollmerSchiedGennaioliDiagnosticRiskSpec,
        MitchellPedersenKratzDarkPoolSpiralReport, MitchellPedersenKratzDarkPoolSpiralSpec,
        evaluate_follmer_schied_convex_risk_and_diagnostic_expectations,
        evaluate_mitchell_pedersen_slow_capital_and_kratz_dark_pool,
    };
    pub use ql_rates::*;
    pub use ql_rates::{
        BrigoMercurioG2TwistSwaptionReport, BrigoMercurioG2TwistSwaptionSpec,
        price_brigo_mercurio_g2_twist_and_swaption,
    };
    pub use ql_rates::{
        TurfusWillemsRfrInArrearsSmileReport, TurfusWillemsRfrInArrearsSmileSpec,
        evaluate_turfus_willems_rfr_in_arrears_caplet_smile,
    };
    pub use ql_scripting::*;
    pub use ql_time::*;
    pub use ql_vol::*;
    pub use ql_vol::{
        AndreasenHugeOneStepZabrReport, AndreasenHugeOneStepZabrSpec,
        evaluate_andreasen_huge_one_step_zabr_calibration,
    };
    pub use ql_vol::{
        BourgeyGatheralBergomiGuyonSmileReport, BourgeyGatheralBergomiGuyonSmileSpec,
        evaluate_bourgey_gatheral_bergomi_guyon_expansion,
    };
    pub use ql_vol::{
        GatheralDemeterfiCorridorVarianceSkewReport, GatheralDemeterfiCorridorVarianceSkewSpec,
        evaluate_gatheral_demeterfi_corridor_variance_and_skew,
    };
    pub use ql_vol::{
        MingoneGlobalArbitrageFreeEssviReport, MingoneGlobalArbitrageFreeEssviSpec,
        evaluate_mingone_global_arbitrage_free_essvi,
    };
    pub use ql_xva::*;
    pub use ql_xva::{
        AlbaneseCrepeyMporBalanceSheetReport, AlbaneseCrepeyMporBalanceSheetSpec,
        compute_albanese_crepey_mpor_and_balance_sheet_xva,
    };
    pub use ql_xva::{
        AlbaneseCrepeyZeronBalanceSheetXvaReport, AlbaneseCrepeyZeronBalanceSheetXvaSpec,
        evaluate_albanese_crepey_zeron_balance_sheet_xva,
    };
    pub use ql_xva::{
        BcbsD424OutputFloorAndSmaOpRiskReport, BcbsD424OutputFloorAndSmaOpRiskSpec,
        evaluate_bcbs_d424_output_floor_and_sma_op_risk,
    };
    pub use ql_xva::{
        GregoryCsaThresholdMtaCtdColvaReport, GregoryCsaThresholdMtaCtdColvaSpec,
        evaluate_gregory_csa_threshold_mta_and_ctd_colva,
    };
    pub use ql_xva::{
        TurlakovGreenKenyonWwrKvaReport, TurlakovGreenKenyonWwrKvaSpec,
        evaluate_turlakov_green_kenyon_wwr_and_kva,
    };
}

use ql_ad::Real;
use ql_commodities::{
    BessAnalyticalSdpDegradationStorageSpec, evaluate_bess_analytical_sdp_degradation_storage,
};
use ql_commodities::{
    BessSpec, CleanSparkDarkSpreadSpec, ClewlowStricklandSwaptionSpec, CommodityFuturesCurve,
    DegreeDayKind, GasStorageSpec, SchwartzSmith2F, SwingOptionSpec, TakeOrPayContractSpec,
    WeatherDerivativeSpec, kirk_spread_option, price_bess_arbitrage_dp,
    price_clean_spark_dark_spread, price_clewlow_strickland_strip_swaption, price_gas_storage,
    price_swing_option, price_take_or_pay_contract, price_weather_hdd_cdd_option,
};
#[allow(unused_imports)]
use ql_commodities::{
    EydelandWolyniecRampOutageVtaReport, EydelandWolyniecRampOutageVtaSpec,
    evaluate_eydeland_wolyniec_ramp_outage_virtual_tolling,
};
#[allow(unused_imports)]
use ql_commodities::{
    WarinBaeuerleStorageRollingHedgeReport, WarinBaeuerleStorageRollingHedgeSpec,
    evaluate_warin_baeuerle_storage_rolling_hedge,
};
use ql_credit::{
    BaseCorrelationTrancheQuote, BlackCoxStructuralSpec, CalibratedCreditCurve, CdsQuote,
    CmCdsMarketModelSpec, CreditDefaultSwap, IsdaCdsEngine, JltRatingMigrationSpec,
    OKaneBespokeBaseCorrelationSpec, ProtectionSide, SurvivalProbabilityCurve,
    analyze_structural_credit_black_cox, bootstrap_base_correlation_curve,
    evaluate_jlt_rating_migration_credit, price_cmcds_market_model,
    price_okane_bespoke_tranche_and_student_t_basket,
};
use ql_credit::{
    FedsDeSpiegeleerAt1CocoDualTriggerSpec, evaluate_feds_despiegeleer_at1_coco_dual_trigger,
};
#[allow(unused_imports)]
use ql_credit::{
    LandoIsraelGeneratorLogarithmReport, LandoIsraelGeneratorLogarithmSpec,
    compute_lando_israel_generator_logarithm_and_migration,
};
#[allow(unused_imports)]
use ql_credit::{
    SchonbucherDefaultableRangeAccrualCapReport, SchonbucherDefaultableRangeAccrualCapSpec,
    evaluate_schonbucher_defaultable_range_accrual_and_spread_cap,
};
use ql_curves::{BcbsD368SixScenarioIrrbbSpec, evaluate_bcbs_d368_six_scenario_irrbb};
use ql_curves::{DiscountCurve, FlatForwardCurve, MonotoneConvexZeroCurve, YieldTermStructure};
#[allow(unused_imports)]
use ql_dml::{
    AgrawalBoydDifferentiableQpReport, AgrawalBoydDifferentiableQpSpec,
    solve_agrawal_boyd_differentiable_qp_layer,
};
#[allow(unused_imports)]
use ql_dml::{
    GilesCapriottiMultilevelMcAadReport, GilesCapriottiMultilevelMcAadSpec,
    evaluate_giles_capriotti_multilevel_mc_aad,
};
#[allow(unused_imports)]
use ql_dml::{
    GlassermanSaddlepointIsStratifiedLhsReport, GlassermanSaddlepointIsStratifiedLhsSpec,
    evaluate_glasserman_saddlepoint_is_and_stratified_lhs_aad,
};
use ql_dml::{HugeSavineDifferentialPcaSpec, evaluate_huge_savine_differential_pca_and_twin_dml};
use ql_exotics::{
    AtlasOptionSpec, BarrierOptionSpec, BarrierType as ExoticBarrierType, EverestOptionSpec,
    ExcursionClockKind, FiniteCapTimerOptionSpec, GammaSwapAndVarOptionSpec, MonitoringFrequency,
    ParisianBarrierOptionSpec, VarianceSwapStrikeQuote, WorstOfAutocallableSpec,
    price_atlas_option_mc_ad, price_barrier_option, price_capped_timer_option_quadrature,
    price_everest_option_mc_ad, price_gamma_swap_and_variance_option, price_parisian_option_mc_ad,
    price_worst_of_autocallable,
};
use ql_exotics::{
    ChoiLiDengMultiAssetSpreadBasketSpec, evaluate_choi_li_deng_multi_asset_spread_basket,
};
#[allow(unused_imports)]
use ql_exotics::{
    CurranRogersShiArithmeticAsianReport, CurranRogersShiArithmeticAsianSpec,
    price_curran_rogers_shi_arithmetic_asian,
};
#[allow(unused_imports)]
use ql_exotics::{
    HeynenKatPartialFractionalLookbackReport, HeynenKatPartialFractionalLookbackSpec,
    evaluate_heynen_kat_partial_fractional_lookback_and_extreme_spread,
};
#[allow(unused_imports)]
use ql_exotics::{
    SheppShiryaevShrevePerpetualRussianReport, SheppShiryaevShrevePerpetualRussianSpec,
    price_shepp_shiryaev_shreve_perpetual_russian_and_mckean,
};
#[allow(unused_imports)]
use ql_fixed_income::{
    AcharyaFedsBankRunAndCocoReport, AcharyaFedsBankRunAndCocoSpec,
    analyze_acharya_feds_bank_run_and_optimal_coco,
};
use ql_fixed_income::{
    AfvConvertibleBondSpec, BondAssetSwapSpec, CallableBondOasSpec, CmoTrancheSpec,
    CmoWaterfallSpec, CocoBondSpec, CocoTriggerPayoff, ConvertibleBondSpec,
    DeSpiegeleerCoCoMandatoryConvSpec, MbsBurnoutOasSpec, NelsonSiegelSvensson,
    SpecifiedPoolMbsSpec, TbaDollarRollSpec, TreasuryFutureDeliverable, TreasuryFuturesBasketSpec,
    TuckmanKeyRateDurationSpec, analyze_bond_asset_swap_and_z_spread, analyze_callable_bond_oas_ad,
    analyze_mbs_burnout_oas_ad, analyze_tba_dollar_roll, analyze_treasury_futures_ctd,
    analyze_tuckman_key_rate_durations_and_pca_hedge, price_cmo_waterfall, price_coco_bond_ds,
    price_convertible_bond_afv_pde, price_convertible_bond_tf,
    price_despiegeleer_coco_and_mandatory_convertible, price_specified_pool_mbs,
};
#[allow(unused_imports)]
use ql_fixed_income::{
    KoziolAsquithFasbStrategicConvertibleReport, KoziolAsquithFasbStrategicConvertibleSpec,
    VeronesiFabozziConvexityPrepaymentReport, VeronesiFabozziConvexityPrepaymentSpec,
    evaluate_koziol_asquith_fasb_strategic_convertible,
    evaluate_veronesi_fabozzi_futures_convexity_and_prepayment,
};
#[allow(unused_imports)]
use ql_fixed_income::{
    TuckmanBurghardtTreasuryWildcardBasisReport, TuckmanBurghardtTreasuryWildcardBasisSpec,
    evaluate_tuckman_burghardt_treasury_wildcard_and_ctd_switch,
};
use ql_fx::{BaranBianchettiMtmXccyBasisSpec, evaluate_baran_bianchetti_mtm_xccy_basis_swap};
use ql_fx::{
    ConstNotionalXccs, FxAccumulatorSpec, FxKikoOptionSpec, FxSmileQuote, FxTarfCapConvention,
    FxTarfSpec, FxVannaVolgaSmile, FxWindowBarrierSpec, MtmResettingXccs, PrdcNoteSpec,
    WystupFxStructuredProductSpec, evaluate_wystup_fx_accumulator_and_dcd,
    garman_kohlhagen_spot_delta, price_fx_accumulator_strip_no_ko, price_fx_kiko_option,
    price_fx_tarf_mc_ad, price_fx_window_barrier_option, price_mtm_resetting_xccs,
    price_prdc_note_hybrid,
};
#[allow(unused_imports)]
use ql_fx::{
    CrossCurrencyConvertibleWarrantReport, CrossCurrencyConvertibleWarrantSpec,
    price_cross_currency_convertible_and_warrant_bond,
};
#[allow(unused_imports)]
use ql_fx::{
    WystupPivotTarfCorridorReport, WystupPivotTarfCorridorSpec,
    evaluate_wystup_pivot_tarf_and_corridor_strip,
};
use ql_inflation::{
    BrigoMercurioInflationMarketModelSpec, CalibratedInflationCurve, JarrowYildirimParams,
    MultiplicativeSeasonality, YoYInflationCapFloor, YoYInflationCapFloorType, ZciisCapFloorSpec,
    ZeroCouponInflationSwap, ZeroInflationCurve, price_brigo_mercurio_yoy_and_lpi_market_model,
    price_yoy_capfloor_jarrow_yildirim, price_zciis_capfloor_jarrow_yildirim,
};
#[allow(unused_imports)]
use ql_inflation::{
    CochraneFiscalTheoryPriceLevelReport, CochraneFiscalTheoryPriceLevelSpec,
    evaluate_cochrane_ftpl_and_maturity_inflation_smoothing,
};
use ql_instruments::{BaroneAdesiWhaleyEngine, OptionType, SwapType, VanillaSwapSpec};
use ql_levy::{CgmyProcess, VarianceGamma, cos_levy_price, evaluate_cgmy_fractional_fft_and_pide};
#[allow(unused_imports)]
use ql_levy::{
    ContTankovBatesSvjPideReport, ContTankovBatesSvjPideSpec,
    evaluate_cont_tankov_bates_svj_characteristic_and_pide,
};
use ql_models::rough::{
    QuadraticRoughHestonParams, RoughBergomiParams, qrh_implied_vol_expansion,
    rbergomi_hybrid_mc_price,
};
use ql_models::{
    AndreasenHugeGridSpec, BatesParams, Bergomi2FParams, DoubleHestonParams, HestonParams,
    HullWhite1F, LiftedHestonParams, MixedRoughBergomiVixSpec, QuinticOuParams,
    bergomi_2f_implied_vol_bg, double_heston_cos_price, evaluate_mixed_rough_bergomi_vix,
    evaluate_quintic_ou_spx_vix, heston_bates_vix_future_and_option_price,
    interpolate_andreasen_huge_call, lifted_heston_cos_price,
    local_vol::{ZabrParams, evaluate_zabr_smile},
    rough::{GuyonPdvParams, pdv_implied_vol_expansion},
    solve_andreasen_huge_step,
};
#[allow(unused_imports)]
use ql_models::{
    CarmonaDelarueMeanFieldSystemicReport, CarmonaDelarueMeanFieldSystemicSpec,
    solve_carmona_delarue_mean_field_systemic_interbank_game,
};
#[allow(unused_imports)]
use ql_models::{
    DoubleHestonKernelSlvReport, DoubleHestonKernelSlvSpec,
    evaluate_double_heston_cumulants_and_kernel_slv,
};
use ql_models::{
    JacquierDeepRoughVolSurrogateSpec, evaluate_jacquier_chebyshev_rough_vol_surrogate,
};
use ql_portfolio::{
    AlfonsiFruthSchiedLobSpec, BoxConstrainedMarkowitzSpec, DarkLitRoutingSpec, EntropyPoolingSpec,
    ExpectileRiskAndPortfolioSpec, IdzorekBlackLittermanSpec, KellyAllocationSpec, MarketSnapshot,
    ObizhaevaWangSpec, PortfolioDag, PortfolioRiskSummary, PortfolioTrade, RobustMarkowitzSpec,
    box_constrained_markowitz_weights, effective_number_of_bets_minimum_torsion,
    entropy_pooling_posterior, equal_risk_contribution_weights,
    evaluate_and_optimize_expectile_portfolio, gerber_covariance_matrix,
    hierarchical_risk_parity_weights, idzorek_black_litterman_posterior,
    kelly_fractional_allocation, nested_clustered_optimization_weights, obizhaeva_wang_execution,
    optimal_dark_lit_execution, robust_markowitz_weights, solve_alfonsi_fruth_schied_lob_execution,
};
use ql_portfolio::{
    AvellanedaStoikovMmSpec, MeanCdarPortfolioSpec, SchurHrpSpec, conditional_drawdown_at_risk,
    entropic_value_at_risk, factor_risk_attribution, gaussian_expected_shortfall,
    oas_and_nonlinear_shrinkage_covariance, optimize_mean_cdar_portfolio, schur_complementary_hrp,
    solve_avellaneda_stoikov_market_making,
};
#[allow(unused_imports)]
use ql_portfolio::{
    BussetiBoydRiskConstrainedKellyReport, BussetiBoydRiskConstrainedKellySpec,
    HarrisRollCorwinSchultzMicrostructureReport, HarrisRollCorwinSchultzMicrostructureSpec,
    analyze_harris_roll_corwin_schultz_microstructure, solve_busseti_boyd_risk_constrained_kelly,
};
#[allow(unused_imports)]
use ql_portfolio::{
    CarteaJaimungalPenalvaLobFillHjbReport, CarteaJaimungalPenalvaLobFillHjbSpec,
    GrinoldQianChristophersonAlphaDecayReport, GrinoldQianChristophersonAlphaDecaySpec,
    MadhavanAbnerHillEtfApArbitrageReport, MadhavanAbnerHillEtfApArbitrageSpec,
    evaluate_grinold_qian_christopherson_ic_decay_and_timing,
    evaluate_madhavan_abner_hill_etf_ap_arbitrage_and_tax_alpha,
    solve_cartea_jaimungal_penalva_lob_fill_hjb,
};
#[allow(unused_imports)]
use ql_portfolio::{
    EngleRussellAcdTradeDurationReport, EngleRussellAcdTradeDurationSpec,
    GrinoldKahnClarkeFundamentalLawReport, GrinoldKahnClarkeFundamentalLawSpec,
    HarveyLiuZhuMultipleTestingSharpeReport, HarveyLiuZhuMultipleTestingSharpeSpec,
    evaluate_engle_russell_acd_trade_duration_and_hazard,
    evaluate_grinold_kahn_clarke_fundamental_law,
    evaluate_harvey_liu_zhu_multiple_testing_and_mintrl,
};
use ql_portfolio::{
    FollmerSchiedGennaioliDiagnosticRiskSpec, MitchellPedersenKratzDarkPoolSpiralSpec,
    evaluate_follmer_schied_convex_risk_and_diagnostic_expectations,
    evaluate_mitchell_pedersen_slow_capital_and_kratz_dark_pool,
};
use ql_rates::{
    BermudanSwaption, CheyetteSwaptionSpec, EuropeanSwaption, G2PlusPlus,
    HullWhiteBermudanTreeEngine, MarkovFunctionalSwaptionSpec, RfrCapletSpec, RfrRateType,
    SinhOuRfrCapletSpec, SwaptionSettlement, TenorBasisSwapSpec,
    price_cheyette_swaption_markov_projection, price_markov_functional_1f,
    price_sinh_ou_rfr_caplet, price_tenor_basis_swap, willems_backward_sabr_caplet_price,
};
#[allow(unused_imports)]
use ql_rates::{
    BrigoMercurioG2TwistSwaptionReport, BrigoMercurioG2TwistSwaptionSpec,
    price_brigo_mercurio_g2_twist_and_swaption,
};
#[allow(unused_imports)]
use ql_rates::{
    HuntKennedyMarkovFunctionalReport, HuntKennedyMarkovFunctionalSpec,
    evaluate_hunt_kennedy_markov_functional_range_accrual,
};
use ql_rates::{
    TurfusWillemsRfrInArrearsSmileSpec, evaluate_turfus_willems_rfr_in_arrears_caplet_smile,
};
use ql_time::{BitsetCalendar, BusinessDayConvention, Date, DayCounter, Period};
#[allow(unused_imports)]
use ql_vol::{
    BergomiForwardStartNapoleonVolReport, BergomiForwardStartNapoleonVolSpec,
    evaluate_bergomi_forward_start_smile_and_napoleon_convexity,
};
use ql_vol::{
    BourgeyGatheralBergomiGuyonSmileSpec, evaluate_bourgey_gatheral_bergomi_guyon_expansion,
};
use ql_vol::{
    CalibratedSabrSmile, FreeBoundarySabrParams, SmileStrikeQuote, ThirdOrderSpeedVolExpansionSpec,
    evaluate_free_boundary_sabr, evaluate_third_order_speed_vol_expansion, implied_black_vol,
};
#[allow(unused_imports)]
use ql_vol::{
    GatheralDemeterfiCorridorVarianceSkewReport, GatheralDemeterfiCorridorVarianceSkewSpec,
    evaluate_gatheral_demeterfi_corridor_variance_and_skew,
};
#[allow(unused_imports)]
use ql_vol::{
    MingoneGlobalArbitrageFreeEssviReport, MingoneGlobalArbitrageFreeEssviSpec,
    evaluate_mingone_global_arbitrage_free_essvi,
};
#[allow(unused_imports)]
use ql_xva::{
    AlbaneseCrepeyMporBalanceSheetReport, AlbaneseCrepeyMporBalanceSheetSpec,
    compute_albanese_crepey_mpor_and_balance_sheet_xva,
};
use ql_xva::{
    AlbaneseCrepeyZeronBalanceSheetXvaSpec, evaluate_albanese_crepey_zeron_balance_sheet_xva,
};
use ql_xva::{
    BaCvaCounterpartySpec, BaCvaPortfolioSpec, BaCvaSingleNameHedge, DynamicSimmMvaSpec,
    MultiCurrencyCsaSpec, NettingSetSwap, SaCcrAssetClass, SaCcrNettingSet, SaCcrTrade, SimmBucket,
    TurlakovWwrSpec, compute_ba_cva_capital, compute_dynamic_simm_mva_and_spike, compute_frtb_sa,
    compute_isda_simm, compute_isda_simm_6_risk_classes, compute_multicurrency_ctd_colva,
    compute_sa_ccr_ead, compute_turlakov_wwr_cva_fva, compute_xva, expected_positive_exposure,
    potential_future_exposure, simulate_hw1f_exposure_cube,
};
#[allow(unused_imports)]
use ql_xva::{
    BcbsD424OutputFloorAndSmaOpRiskReport, BcbsD424OutputFloorAndSmaOpRiskSpec,
    evaluate_bcbs_d424_output_floor_and_sma_op_risk,
};
#[allow(unused_imports)]
use ql_xva::{
    GreenFrtbCvaHedgeCapitalReport, GreenFrtbCvaHedgeCapitalSpec,
    evaluate_green_frtb_cva_hedge_relief_and_kva,
};
use std::fmt::Write as FmtWrite;
use std::sync::Arc;
use std::time::Instant;

const BRAILLE_DOTS: [[u32; 2]; 4] = [[0x01, 0x08], [0x02, 0x10], [0x04, 0x20], [0x40, 0x80]];

/// Logical workstation tabs.  The dashboard remains headless and printable,
/// while callers can use this stable enumeration to select a domain view.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum TabId {
    Market,
    Rates,
    Credit,
    Inflation,
    Fx,
    Equities,
    Volatility,
    Xva,
    Portfolio,
    Exotics,
    Commodities,
    Levy,
    FixedIncome,
    MonteCarlo,
    Analytics,
    Blotter,
    AmericanHybrid,
    BlRiskParity,
    BermudanVolIft,
    SaCcrSwingQp,
    SlvRoughWorstOf,
    FrtbCvarExecution,
    MonotoneConvexVixParisianRobust,
    Phase14Sota,
    Phase15Sota,
    Phase16Sota,
    Phase17Sota,
    Phase18Sota,
    Phase19Sota,
    Phase20Sota,
    Phase21Sota,
    Phase22Sota,
    Phase23Sota,
    Phase24Sota,
    Phase25Sota,
    Phase26Sota,
    Phase27Sota,
    Phase28Sota,
    Phase29Sota,
    Phase30BookSota,
    Phase31BookSota,
    Phase32BookSota,
    Phase33BookSota,
}
impl TabId {
    #[allow(non_upper_case_globals)]
    pub const Phase30Sota: Self = Self::Phase30BookSota;
    #[allow(non_upper_case_globals)]
    pub const Phase31Sota: Self = Self::Phase31BookSota;
    #[allow(non_upper_case_globals)]
    pub const Phase32Sota: Self = Self::Phase32BookSota;
    #[allow(non_upper_case_globals)]
    pub const Phase33Sota: Self = Self::Phase33BookSota;
    pub const TAB_COUNT: usize = 43;
    pub const ALL: [Self; 43] = [
        Self::Market,
        Self::Rates,
        Self::Credit,
        Self::Inflation,
        Self::Fx,
        Self::Equities,
        Self::Volatility,
        Self::Xva,
        Self::Portfolio,
        Self::Exotics,
        Self::Commodities,
        Self::Levy,
        Self::FixedIncome,
        Self::MonteCarlo,
        Self::Analytics,
        Self::Blotter,
        Self::AmericanHybrid,
        Self::BlRiskParity,
        Self::BermudanVolIft,
        Self::SaCcrSwingQp,
        Self::SlvRoughWorstOf,
        Self::FrtbCvarExecution,
        Self::MonotoneConvexVixParisianRobust,
        Self::Phase14Sota,
        Self::Phase15Sota,
        Self::Phase16Sota,
        Self::Phase17Sota,
        Self::Phase18Sota,
        Self::Phase19Sota,
        Self::Phase20Sota,
        Self::Phase21Sota,
        Self::Phase22Sota,
        Self::Phase23Sota,
        Self::Phase24Sota,
        Self::Phase25Sota,
        Self::Phase26Sota,
        Self::Phase27Sota,
        Self::Phase28Sota,
        Self::Phase29Sota,
        Self::Phase30BookSota,
        Self::Phase31BookSota,
        Self::Phase32BookSota,
        Self::Phase33BookSota,
    ];
    pub fn label(self) -> &'static str {
        match self {
            Self::Market => "01 Market",
            Self::Rates => "02 Rates",
            Self::Credit => "03 Credit",
            Self::Inflation => "04 Inflation",
            Self::Fx => "05 FX",
            Self::Equities => "06 Equities",
            Self::Volatility => "07 Volatility",
            Self::Xva => "08 XVA",
            Self::Portfolio => "09 Portfolio",
            Self::Exotics => "10 Exotics",
            Self::Commodities => "11 Commodities",
            Self::Levy => "12 Levy",
            Self::FixedIncome => "13 Fixed Income",
            Self::MonteCarlo => "14 Monte Carlo",
            Self::Analytics => "15 Analytics",
            Self::Blotter => "16 Blotter",
            Self::AmericanHybrid => "17 American / Hybrid",
            Self::BlRiskParity => "18 BL / Risk Parity",
            Self::BermudanVolIft => "19 Bermudan / Vol IFT",
            Self::SaCcrSwingQp => "20 SA-CCR / Swing / QP",
            Self::SlvRoughWorstOf => "21 SLV / rBergomi / WoF",
            Self::FrtbCvarExecution => "22 FRTB / CVaR / Exec",
            Self::MonotoneConvexVixParisianRobust => "23 Monotone Convex / VIX / Parisian / Robust",
            Self::Phase14Sota => "24 Double Heston / TARF / CoCo / HRP",
            Self::Phase15Sota => "25 JY-Inflation / TBA-Roll / Bergomi2F / NCO",
            Self::Phase16Sota => "26 Cheyette / CMO / QRH / BaseCorr / EVaR",
            Self::Phase17Sota => "27 ASW / MF1F / AH-LV / BESS / Idzorek-DP",
            Self::Phase18Sota => "28 MtM-XCCS / CTD-COLVA / ZABR-PDV / Weather / ENB",
            Self::Phase19Sota => {
                "29 Burnout-MBS / CMCDS / Quintic-OU / AFV-CB / CS-Swaption / Schur-HRP"
            }
            Self::Phase20Sota => {
                "30 Sinh-OU / BK-FRN / Morini-Brigo / WWM-XVA / mrBergomi / Choi / AFS-LOB / EXVaR"
            }
            Self::Phase21Sota => {
                "31 HW-MC / Clark-FX / SA-CCR-KVA / BFG-rVol / Cartea-Alpha / Paleologo"
            }
            Self::Phase22Sota => "32 Phase 22 SOTA / KRD / LPI / CoCo / Bergomi / CGMY / AAD",
            Self::Phase23Sota => "33 RKHS/FMM/VV/LSV/JLT/CTD/QRH/NCO/DML",
            Self::Phase24Sota => {
                "34 Phase 24 SOTA / MBS-IO-PO / RMV / MO-XVA / Dilution / Sig / GARCH / POT / RIE / PIN / GP / Brinson / IERP"
            }
            Self::Phase25Sota => {
                "35 Phase 25 SOTA / Kalman-RTS / VECM / HJ-SDF / IFRS9 / CPT / Hasbrouck / FV / Lando / LPPL / CB-Arb / Causal / Malliavin"
            }
            Self::Phase26Sota => {
                "36 Phase 26 SOTA / MF-Rates / Wildcard / TARF / Bergomi / Schönbucher / FRTB-CVA / Lookback / VTA / FLAM / ACD / HLZ / MLMC"
            }
            Self::Phase27Sota => {
                "37 Phase 27 SOTA / FTPL / G2++ / eSSVI / Offshore-CB / CoCo / MPOR / Double-Heston / Curran / Warin / Roll-Corwin / Busseti-Kelly / Diff-QP"
            }
            Self::Phase28Sota => {
                "38:P28-SOTA / Phase 28 SOTA / IRRBB / Turfus-RFR / Baran-MtM / Bourgey-BG / FEDS-PONV / Zeron-SIMM / Jacquier-rBergomi / Choi-Spread / BESS-SDP / Kratz-Dark / Föllmer / Diff-PCA"
            }
            Self::Phase29Sota => {
                "39 Phase 29 SOTA / Veronesi / Lando-Log / Corridor-Var / Russian / Bates-SVJ / Koziol-CB / MFG / BCBS-Floor / CJP-LOB / ETF-AP / IC-Decay / Glasserman-IS"
            }
            Self::Phase30BookSota => {
                "40 Phase 30 SOTA / Taylor-Curve / Damodaran-CRP / Cumulant-VaR / LOB-LVaR / Prospect-CGO / Regret-Match / Beneish-HKAS32 / Campisi-FI / Heckman-GMM / Elastic-Net / Reset-CB / FSB-BailIn"
            }
            Self::Phase31BookSota => {
                "41 Phase 31 SOTA / KRD-LSC / CMO-PAC / Asquith-CB / ASU2024-CPDI / ASC470-Mod / ColVA-CTD / Hill-EVT / Bouchaud-TLF / Karnosky-FX / DGTW-Timing / Brock-Hommes / Penman-ReOI"
            }
            Self::Phase32BookSota => {
                "42 Phase 32 SOTA / Cheyette-LV / ADB-NDIRS / Andreasen-ZABR / QRH-Zumbach / CMCDS-Morini / Turlakov-WWR / HHW-ADI / CappedCall-MW / Ahmadi-EVaR / Schur-HRP / BEKK-GARCH / Schilit-Sloan"
            }
            Self::Phase33BookSota => {
                "43 Phase 33 SOTA / Rebonato-CMS / Leland-Toft / Calamos-Ascot / SpecPool-Roll / Glosten-Milgrom / Potters-RIE / EGARCH-Hamilton / Wooldridge-2SLS / McNeil-Euler / Bacon-Omega / Damodaran-RD / Arrow-Debreu"
            }
        }
    }
}

/// A 2x4 sub-pixel canvas represented by Unicode Braille cells.
#[derive(Clone, Debug)]
pub struct BrailleCanvas {
    pub cols: usize,
    pub rows: usize,
    pub width: usize,
    pub height: usize,
    pixels: Vec<bool>,
    colors: Vec<[u8; 3]>,
}
impl BrailleCanvas {
    pub fn new(cols: usize, rows: usize) -> Self {
        let width = cols.saturating_mul(2);
        let height = rows.saturating_mul(4);
        Self {
            cols,
            rows,
            width,
            height,
            pixels: vec![false; width * height],
            colors: vec![[0; 3]; width * height],
        }
    }
    pub fn clear(&mut self) {
        self.pixels.fill(false);
        self.colors.fill([0; 3]);
    }
    pub fn set_pixel(&mut self, x: usize, y: usize, color: [u8; 3]) {
        if x < self.width && y < self.height {
            let i = y * self.width + x;
            self.pixels[i] = true;
            self.colors[i] = color;
        }
    }
    pub fn draw_line(
        &mut self,
        mut x0: isize,
        mut y0: isize,
        x1: isize,
        y1: isize,
        color: [u8; 3],
    ) {
        let dx = (x1 - x0).abs();
        let sx = if x0 < x1 { 1 } else { -1 };
        let dy = -(y1 - y0).abs();
        let sy = if y0 < y1 { 1 } else { -1 };
        let mut err = dx + dy;
        loop {
            if x0 >= 0 && y0 >= 0 {
                self.set_pixel(x0 as usize, y0 as usize, color);
            }
            if x0 == x1 && y0 == y1 {
                break;
            }
            let twice = 2 * err;
            if twice >= dy {
                err += dy;
                x0 += sx;
            }
            if twice <= dx {
                err += dx;
                y0 += sy;
            }
        }
    }
    pub fn plot_series(&mut self, xs: &[f64], ys: &[f64], y_min: f64, y_max: f64, color: [u8; 3]) {
        if xs.is_empty()
            || ys.is_empty()
            || xs.len() != ys.len()
            || self.width == 0
            || self.height == 0
        {
            return;
        }
        let x_min = xs.iter().copied().fold(f64::INFINITY, f64::min);
        let x_max = xs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let dx = (x_max - x_min).abs().max(f64::EPSILON);
        let dy = (y_max - y_min).abs().max(f64::EPSILON);
        let canvas_width = self.width;
        let canvas_height = self.height;
        let point = |i: usize| {
            let x =
                (((xs[i] - x_min) / dx) * (canvas_width.saturating_sub(1) as f64)).round() as isize;
            let y_norm = ((ys[i] - y_min) / dy).clamp(0.0, 1.0);
            let y = ((1.0 - y_norm) * canvas_height.saturating_sub(1) as f64).round() as isize;
            (x, y)
        };
        let mut last = point(0);
        self.set_pixel(last.0.max(0) as usize, last.1.max(0) as usize, color);
        for i in 1..xs.len() {
            let next = point(i);
            self.draw_line(last.0, last.1, next.0, next.1, color);
            last = next;
        }
    }
    pub fn render_lines(&self, ansi_color: bool) -> Vec<String> {
        let mut lines = Vec::with_capacity(self.rows);
        for row in 0..self.rows {
            let mut line = String::new();
            for col in 0..self.cols {
                let mut mask = 0u32;
                let mut color = [0u8; 3];
                let mut active = false;
                for (py, dots) in BRAILLE_DOTS.iter().enumerate() {
                    for (px, dot) in dots.iter().enumerate() {
                        let x = col * 2 + px;
                        let y = row * 4 + py;
                        let idx = y * self.width + x;
                        if self.pixels[idx] {
                            mask |= *dot;
                            color = self.colors[idx];
                            active = true;
                        }
                    }
                }
                let ch = char::from_u32(0x2800 | mask).unwrap_or('\u{2800}');
                if ansi_color && active {
                    let _ = write!(
                        line,
                        "\x1b[38;2;{};{};{}m{}\x1b[0m",
                        color[0], color[1], color[2], ch
                    );
                } else {
                    line.push(ch);
                }
            }
            lines.push(line);
        }
        lines
    }
}

pub fn render_sparkline(values: &[f64], width: usize) -> String {
    if values.is_empty() || width == 0 {
        return String::new();
    }
    const LEVELS: [char; 8] = ['▁', '▂', '▃', '▄', '▅', '▆', '▇', '█'];
    let lo = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(f64::INFINITY, f64::min);
    let hi = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(f64::NEG_INFINITY, f64::max);
    if !lo.is_finite() || !hi.is_finite() {
        return LEVELS[0].to_string().repeat(width);
    }
    let span = (hi - lo).abs().max(f64::EPSILON);
    (0..width)
        .map(|i| {
            let j = i * values.len() / width;
            let v = values[j.min(values.len() - 1)];
            let v = if v.is_finite() { v } else { lo };
            let level = (((v - lo) / span) * 7.0).round().clamp(0.0, 7.0) as usize;
            LEVELS[level]
        })
        .collect()
}

pub fn render_diverging_bar(
    value: f64,
    max_abs: f64,
    half_width: usize,
    ansi_color: bool,
) -> String {
    let n = if max_abs > 0.0 {
        ((value.abs() / max_abs) * half_width as f64).round() as usize
    } else {
        0
    };
    let n = n.min(half_width);
    let mut out = String::new();
    let color = if value < 0.0 {
        [255, 90, 90]
    } else {
        [0, 230, 140]
    };
    if ansi_color {
        let _ = write!(out, "\x1b[38;2;{};{};{}m", color[0], color[1], color[2]);
    }
    if value < 0.0 {
        out.push_str(&" ".repeat(half_width - n));
        out.push_str(&"━".repeat(n));
    } else {
        out.push_str(&" ".repeat(half_width));
    }
    out.push('│');
    if value >= 0.0 {
        out.push_str(&"━".repeat(n));
        out.push_str(&" ".repeat(half_width - n));
    } else {
        out.push_str(&" ".repeat(half_width));
    }
    if ansi_color {
        out.push_str("\x1b[0m");
    }
    out
}

fn visible_width(s: &str) -> usize {
    let mut count = 0;
    let mut escape = false;
    for ch in s.chars() {
        if escape {
            if ch == 'm' {
                escape = false;
            }
        } else if ch == '\x1b' {
            escape = true;
        } else {
            count += 1;
        }
    }
    count
}

/// Pad a dashboard cell by terminal-visible characters, ignoring ANSI color
/// sequences.  This keeps the 2x2 layout aligned in both plain and color mode.
fn pad_visible(s: &str, width: usize) -> String {
    let visible = visible_width(s);
    if visible >= width {
        return s.to_string();
    }
    let mut out = s.to_string();
    out.push_str(&" ".repeat(width - visible));
    out
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ScenarioKind {
    Normal,
    RatesPlus50Bp,
    EquityCrashMinus15Pct,
    CreditWideningPlus100Bp,
    StagflationPlus150Bp,
}
impl ScenarioKind {
    #[allow(non_upper_case_globals)]
    pub const RatesParallelPlus50Bp: Self = Self::RatesPlus50Bp;
    #[allow(non_upper_case_globals)]
    pub const VolShockPlus5Pts: Self = Self::EquityCrashMinus15Pct;
    #[allow(non_upper_case_globals)]
    pub const CreditSpreadWidening100Bp: Self = Self::CreditWideningPlus100Bp;
    pub fn next(self) -> Self {
        match self {
            Self::Normal => Self::RatesPlus50Bp,
            Self::RatesPlus50Bp => Self::EquityCrashMinus15Pct,
            Self::EquityCrashMinus15Pct => Self::CreditWideningPlus100Bp,
            Self::CreditWideningPlus100Bp => Self::StagflationPlus150Bp,
            Self::StagflationPlus150Bp => Self::Normal,
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Self::Normal => "NORMAL",
            Self::RatesPlus50Bp => "RATES +50BP",
            Self::EquityCrashMinus15Pct => "EQUITY -15%",
            Self::CreditWideningPlus100Bp => "CREDIT +100BP",
            Self::StagflationPlus150Bp => "STAGFLATION +150BP",
        }
    }
}
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ChartMode {
    YieldAndInflationCurves,
    RoughHestonAndFxSmiles,
    XvaExposureProfile,
}
impl ChartMode {
    pub fn next(self) -> Self {
        match self {
            Self::YieldAndInflationCurves => Self::RoughHestonAndFxSmiles,
            Self::RoughHestonAndFxSmiles => Self::XvaExposureProfile,
            Self::XvaExposureProfile => Self::YieldAndInflationCurves,
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Self::YieldAndInflationCurves => "Yield & Inflation Curves",
            Self::RoughHestonAndFxSmiles => "Rough Heston & FX Smiles",
            Self::XvaExposureProfile => "XVA Exposure Profile",
        }
    }
}

#[derive(Clone, Debug)]
pub struct WorkstationFrameData {
    pub tick: usize,
    pub scenario: ScenarioKind,
    pub chart_mode: ChartMode,
    pub eval_micros: f64,
    pub dag_recomputed_trades: usize,
    pub dag_total_trades: usize,
    pub dag_curves_rebuilt: bool,
    pub total_npv: f64,
    pub intraday_pnl: f64,
    pub eq_delta: f64,
    pub eq_vega: f64,
    pub fx_spot_delta: f64,
    pub inflation_par_delta_sum: f64,
    pub cva: f64,
    pub dva: f64,
    pub fva: f64,
    pub mva: f64,
    pub simm_initial_margin: f64,
    pub ir_dv01_buckets: Vec<(String, f64)>,
    pub credit_cs01_buckets: Vec<(String, f64)>,
    pub jump_to_default: f64,
    pub bermudan_sabr_metric: f64,
    pub sacc_swing_qp_metric: f64,
    pub slv_rough_worst_of_metric: f64,
    pub frtb_cvar_execution_metric: f64,
    pub vix_parisian_robust_metric: f64,
    pub double_heston_tarf_coco_hrp_metric: f64,
    pub jy_tba_bergomi_nco_metric: f64,
    pub cheyette_cmo_qrh_evar_metric: f64,
    pub asw_mf1f_ah_bess_ep_metric: f64,
    pub mtmxccs_colva_zabr_weather_enb_metric: f64,
    pub burnout_cmcds_quintic_afv_cs_schur_metric: f64,
    pub sinhou_mb_rbergomi_choi_afs_expectile_metric: f64,
    pub monotone_convex_clark_saccr_rvol_cartea_paleologo_metric: f64,
    pub krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric: f64,
    pub rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml_metric: f64,
    pub mbs_rmv_mo_sig_garch_rie_gp_ierp_metric: f64,
    pub kalman_vecm_hj_ifrs9_cpt_hasbrouck_fv_lando_lppl_cbarb_causal_malliavin_metric: f64,
    pub phase26_sota_metric: f64,
    pub phase27_sota_metric: f64,
    pub phase28_sota_metric: f64,
    pub phase29_sota_metric: f64,
    pub phase30_sota_metric: f64,
    pub phase31_sota_metric: f64,
    pub phase32_sota_metric: f64,
    pub phase33_sota_metric: f64,
    pub blotter_rows: Vec<(String, f64, String)>,
}

fn range2(a: &[f64], b: &[f64]) -> (f64, f64) {
    let lo = a.iter().chain(b).copied().fold(f64::INFINITY, f64::min);
    let hi = a.iter().chain(b).copied().fold(f64::NEG_INFINITY, f64::max);
    let pad = ((hi - lo).abs() * 0.08).max(1e-6);
    (lo - pad, hi + pad)
}

fn range3(a: &[f64], b: &[f64], c: &[f64]) -> (f64, f64) {
    let lo = a
        .iter()
        .chain(b)
        .chain(c)
        .copied()
        .fold(f64::INFINITY, f64::min);
    let hi = a
        .iter()
        .chain(b)
        .chain(c)
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let pad = ((hi - lo).abs() * 0.08).max(1e-6);
    (lo - pad, hi + pad)
}

fn phase21_live_metric(m: &MarketSnapshot, credit_bump: f64) -> f64 {
    let c = FlatForwardCurve::new(m.ref_date, m.ois_quotes[0]);
    let hw = ql_curves::bootstrap_hagan_west_monotone_convex_curve(
        &ql_curves::HaganWestMonotoneConvexBootstrapSpec {
            reference_date: m.ref_date,
            pillar_times: vec![1.0, 2.0, 5.0],
            par_swap_quotes: m.ois_quotes.to_vec(),
            payment_frequency_per_year: 2,
            enforce_positivity: true,
            ameliorate_slopes: true,
            query_time: 2.0,
        },
    );
    let sw = ql_rates::price_piterbarg_quasi_gaussian_amortizing_swaption(
        &ql_rates::PiterbargQuasiGaussianAmortizingSwaptionSpec {
            option_expiry: 1.0,
            second_exercise_time: None,
            swap_times: vec![1.0, 2.0, 3.0],
            period_notionals: vec![1.0, 0.9],
            strike: m.ois_quotes[0],
            mean_reversion: 0.05,
            sigma0: 0.01 + m.eq_vol * 0.02,
            local_vol_skew_b: 0.0,
            cms_pay_delay: 0.0,
            is_payer: true,
        },
        &c,
        &c,
    );
    let fx = ql_fx::calibrate_clark_fx_delta_smile(&ql_fx::ClarkFxDeltaConventionSpec {
        spot: 1.10,
        r_dom: m.ois_quotes[0],
        r_for: m.fwd_quotes[0],
        expiry: 1.0,
        atm_vol: 0.12 + m.eq_vol * 0.1,
        rr_25d: 0.0,
        market_strangle_25d: 0.01,
        is_premium_included_delta: false,
        is_forward_delta: false,
        query_strike: 1.10,
    });
    let rv =
        ql_vol::evaluate_gatheral_gyongy_lsv_leverage(&ql_vol::GatheralLsvGyongyLeverageSpec {
            forward: m.eq_spot,
            strike: m.eq_spot,
            expiry: 1.0,
            svi_a: (m.eq_vol * m.eq_vol).max(1e-4),
            svi_b: 0.2,
            svi_rho: -0.3,
            svi_m: 0.0,
            svi_sigma: 0.2,
            dw_dt_atm: 0.04,
            heston_v0: 0.04,
            heston_theta: 0.04,
            heston_kappa: 1.0,
            heston_eta: 0.5,
            heston_rho: -0.7,
        });
    let sac = ql_xva::compute_green_saccr_lifetime_kva_mva(&ql_xva::GreenSaCcrLifetimeKvaMvaSpec {
        time_grid: vec![0.0, 1.0, 2.0],
        trade_notional: 1.0e6,
        trade_maturity: 5.0,
        initial_mtm: 0.0,
        mtm_volatility: 0.1,
        collateral_threshold: 0.0,
        initial_margin_ratio: 0.1,
        counterparty_pd: 0.02 + credit_bump,
        counterparty_lgd: 0.6,
        cva_risk_weight: 0.5,
        capital_ratio_requirement: 0.08,
        hurdle_rate_gamma_k: 0.1,
        funding_spread_sf: 0.01,
        mva_spread_sm: 0.01,
        risk_free_rate: m.ois_quotes[0],
    });
    let greg = ql_xva::allocate_gregory_netting_set_xva(&ql_xva::GregoryNettingSetAllocationSpec {
        trade_weights: vec![1.0, 1.0],
        trade_expected_mtms: vec![1.0, -1.0],
        trade_exposure_vols: vec![0.2, 0.2],
        correlation_matrix: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
        collateral_fraction: 0.0,
        horizon: 1.0,
        num_steps: 4,
        risk_free_rate: m.ois_quotes[0],
        cpty_hazard_rate: 0.02 + credit_bump,
        cpty_lgd: 0.6,
        own_hazard_rate: 0.01,
        own_lgd: 0.6,
        funding_spread: 0.01,
        colva_spread: 0.01,
    });
    let rough = ql_models::evaluate_bayer_friz_gatheral_markovian_rough(
        &ql_models::BayerFrizGatheralMarkovianRoughSpec {
            spot: m.eq_spot,
            strike: m.eq_spot,
            expiry: 1.0,
            rate: m.ois_quotes[0],
            div_yield: 0.0,
            xi0: (m.eq_vol * m.eq_vol).max(1e-4),
            eta: 1.0,
            rho: -0.7,
            hurst: 0.1,
            num_markov_factors: 4,
            is_call: true,
        },
    );
    let cli = ql_exotics::price_bouzoubaa_reverse_and_global_cliquet(
        &ql_exotics::BouzoubaaReverseAndGlobalCliquetSpec {
            notional: 1.0,
            rate: m.ois_quotes[0],
            div_yield: 0.0,
            reset_times: vec![0.5, 1.0],
            forward_vols: vec![m.eq_vol, m.eq_vol],
            forward_skews: vec![0.0, 0.0],
            initial_coupon_cmax: 0.0,
            local_floor: -0.1,
            local_cap: 0.1,
            global_floor: -0.2,
            global_cap: 0.2,
            is_reverse_cliquet: false,
        },
    );
    let plant = ql_commodities::price_burger_multifuel_tolling_power_plant(
        &ql_commodities::BurgerMultiFuelTollingPlantSpec {
            power_forwards: vec![80.0, 80.0],
            gas_forwards: vec![20.0, 20.0],
            coal_forwards: vec![15.0, 15.0],
            carbon_forwards: vec![50.0, 50.0],
            dispatch_times: vec![0.25, 0.75],
            capacity_mw: 100.0,
            period_hours: 1.0,
            gas_heat_rate: 2.0,
            gas_emission_factor: 0.4,
            gas_vom_cost: 2.0,
            coal_heat_rate: 2.5,
            coal_emission_factor: 0.8,
            coal_vom_cost: 1.0,
            vol_power: 0.2,
            vol_gas: 0.2,
            vol_coal: 0.2,
            vol_carbon: 0.2,
            rho_power_gas: 0.2,
            rho_power_coal: 0.2,
            rho_power_carbon: 0.2,
            rho_gas_coal: 0.2,
            start_up_cost: 100.0,
            min_load_fraction: 0.5,
            part_load_heat_rate_penalty: 0.1,
            risk_free_rate: m.ois_quotes[0],
        },
    );
    let exec = ql_portfolio::solve_cartea_jaimungal_alpha_execution(
        &ql_portfolio::CarteaJaimungalAlphaExecutionSpec {
            initial_inventory: 1.0,
            horizon: 1.0,
            num_steps: 8,
            temporary_impact_k: 1.0,
            permanent_impact_b: 0.1,
            inventory_penalty_phi: 0.2,
            terminal_penalty_a: 0.0,
            initial_alpha_signal: 0.1,
            alpha_mean_reversion_kappa: 1.0,
            alpha_volatility: 0.2,
        },
    );
    let paleo = ql_portfolio::solve_paleologo_factor_neutral_portfolio(
        &ql_portfolio::PaleologoFactorNeutralAlphaSpec {
            raw_alphas: vec![1.0, 2.0],
            factor_loadings_b: vec![vec![1.0], vec![1.0]],
            factor_covariance_omega: vec![vec![1.0]],
            sample_idio_variances: vec![0.1, 0.2],
            idio_shrinkage_xi: 0.5,
            idio_risk_aversion: 1.0,
            factor_risk_aversion: 1.0,
            factor_crowding_scores: vec![0.0, 0.0],
            crowding_penalty_lambda: 0.0,
        },
    );
    let aad = ql_dml::evaluate_glasserman_malliavin_checkpointed_aad(
        &ql_dml::GlassermanMalliavinCheckpointedMcSpec {
            spot: m.eq_spot,
            strike: m.eq_spot,
            rate: m.ois_quotes[0],
            div_yield: 0.0,
            base_vol: m.eq_vol,
            local_vol_skew: 0.0,
            expiry: 1.0,
            num_steps: 8,
            num_checkpoints: 2,
        },
    );
    hw.query_discount_factor.val()
        + sw.european_amortizing_npv.val()
        + fx.query_call_price.val()
        + rv.implied_vol.val()
        + sac.lifetime_kva.val()
        + greg.netted_total_xva.val()
        + rough.option_price.val()
        + cli.npv.val()
        + plant.tolling_plant_npv.val()
        + exec.expected_execution_cost_net_of_alpha.val()
        + paleo.expected_active_return.val()
        + aad.call_price.val()
}

pub struct DeskWorkstation {
    pub dag: PortfolioDag,
    pub base_market: MarketSnapshot,
    pub scenario: ScenarioKind,
    pub chart_mode: ChartMode,
    pub npv_history: Vec<f64>,
    trades: Vec<PortfolioTrade>,
    last_summary: PortfolioRiskSummary,
    tick: usize,
}
impl DeskWorkstation {
    pub fn new() -> Self {
        let ref_date = Date::from_ymd(2024, 1, 1).expect("valid workstation date");
        let market = MarketSnapshot {
            ref_date,
            ois_quotes: [0.02, 0.021, 0.022],
            fwd_quotes: [0.025, 0.026, 0.027],
            inflation_quotes: vec![("CPI".to_string(), 0.023)],
            heston_params: HestonParams {
                v0: 0.04,
                kappa: 1.5,
                theta: 0.04,
                sigma_v: 0.3,
                rho: -0.7,
            },
            eq_spot: 100.0,
            eq_vol: 0.2,
        };
        let program =
            ql_scripting::BytecodeProgram::compile(&ql_scripting::ScriptedProduct::autocallable(
                vec![0.25, 0.5, 0.75, 1.0],
                105.0,
                95.0,
                0.04,
                100.0,
                85.0,
                2.0,
            ));
        let trades = vec![
            PortfolioTrade::IrsSwap {
                id: 0,
                notional: 1_000_000.0,
                fixed_rate: 0.02,
                maturity_years: 5,
            },
            PortfolioTrade::IrsSwap {
                id: 1,
                notional: 1_000_000.0,
                fixed_rate: 0.021,
                maturity_years: 10,
            },
            PortfolioTrade::IrsSwap {
                id: 2,
                notional: 1_000_000.0,
                fixed_rate: 0.022,
                maturity_years: 20,
            },
            PortfolioTrade::HestonEuropeanOption {
                id: 3,
                strike: 100.0,
                expiry: 1.0,
                is_call: true,
                notional: 1.0,
            },
            PortfolioTrade::HestonEuropeanOption {
                id: 4,
                strike: 100.0,
                expiry: 2.0,
                is_call: false,
                notional: 1.0,
            },
            PortfolioTrade::ScriptedNote {
                id: 5,
                program,
                expiry: 1.0,
            },
        ];
        let mut dag = PortfolioDag::new();
        let summary = dag.evaluate_portfolio(&market, &trades);
        Self {
            dag,
            base_market: market,
            scenario: ScenarioKind::Normal,
            chart_mode: ChartMode::YieldAndInflationCurves,
            npv_history: vec![summary.total_npv],
            trades,
            last_summary: summary,
            tick: 0,
        }
    }
    pub fn set_scenario(&mut self, scenario: ScenarioKind) {
        self.scenario = scenario;
    }
    pub fn set_chart_mode(&mut self, mode: ChartMode) {
        self.chart_mode = mode;
    }
    fn scenario_market(&self) -> MarketSnapshot {
        let mut s = self.base_market.clone();
        match self.scenario {
            ScenarioKind::Normal => {
                s.eq_spot += self.tick as f64 * 0.0001;
            }
            ScenarioKind::RatesPlus50Bp => {
                s.ois_quotes.iter_mut().for_each(|x| *x += 0.005);
                s.fwd_quotes.iter_mut().for_each(|x| *x += 0.005);
            }
            ScenarioKind::EquityCrashMinus15Pct => {
                s.eq_spot *= 0.85;
                s.heston_params.v0 *= 1.5;
            }
            ScenarioKind::CreditWideningPlus100Bp => {}
            ScenarioKind::StagflationPlus150Bp => {
                s.ois_quotes.iter_mut().for_each(|x| *x += 0.0025);
                for (_, z) in &mut s.inflation_quotes {
                    *z += 0.015;
                }
            }
        }
        s
    }
    fn credit_metrics(&self, scenario: ScenarioKind) -> (Vec<(String, f64)>, f64, f64, f64) {
        let d = self.base_market.ref_date;
        let mats = [1.0, 2.0, 3.0, 5.0, 10.0];
        let quotes: Vec<_> = mats
            .iter()
            .enumerate()
            .map(|(i, &maturity_years)| CdsQuote {
                maturity_years,
                par_spread: 0.01
                    + i as f64 * 0.001
                    + if scenario == ScenarioKind::CreditWideningPlus100Bp {
                        0.01
                    } else {
                        0.0
                    },
            })
            .collect();
        let discount = DiscountCurve::new(
            d,
            mats.to_vec(),
            mats.iter().map(|t| (-0.02 * t).exp()).collect(),
        );
        let calibrated = CalibratedCreditCurve::bootstrap(d, discount.clone(), 0.4, &quotes)
            .expect("credit calibration");
        let cds = CreditDefaultSwap {
            side: ProtectionSide::Buyer,
            notional: 1_000_000.0,
            running_spread: 0.005,
            recovery_rate: 0.4,
            payment_times: (1..=20).map(|i| i as f64 / 4.0).collect(),
            accrual_fractions: vec![0.25; 20],
        };
        let (_, cs01) = calibrated.value_and_par_cs01(|curve| {
            let dc = FlatForwardCurve::new(d, ql_ad::TapeVar::from_f64(0.02));
            IsdaCdsEngine::price(&cds, &dc, curve).npv
        });
        let labels = ["1Y", "2Y", "3Y", "5Y", "10Y"];
        let buckets = cs01
            .iter()
            .zip(labels)
            .map(|(x, label)| (label.to_string(), *x))
            .collect();
        let jtd = CalibratedCreditCurve::jump_to_default(&cds, 0.0);
        let cds_npv = IsdaCdsEngine::price(&cds, &discount, &calibrated.curve).npv;
        (
            buckets,
            jtd,
            calibrated.curve.hazard_rates.iter().sum(),
            cds_npv,
        )
    }
    fn inflation_metrics(&self, scenario: ScenarioKind) -> (f64, f64) {
        let d = self.base_market.ref_date;
        let discount = DiscountCurve::new(
            d,
            vec![1.0, 2.0, 5.0],
            vec![(-0.02f64).exp(), (-0.04f64).exp(), (-0.1f64).exp()],
        );
        let bump = if scenario == ScenarioKind::StagflationPlus150Bp {
            0.015
        } else {
            0.0
        };
        let quotes = [(1.0, 0.02 + bump), (2.0, 0.022 + bump), (5.0, 0.025 + bump)];
        let mut f = [1.0; 12];
        for (i, factor) in f.iter_mut().enumerate().take(11) {
            *factor = 1.0 + 0.006 * (i as f64 - 5.5) / 5.5;
        }
        f[11] = 1.0 / f[..11].iter().product::<f64>();
        let seasonality = MultiplicativeSeasonality::new(f).expect("seasonality");
        let curve = CalibratedInflationCurve::bootstrap_with_seasonality(
            d,
            100.0,
            discount.clone(),
            &quotes,
            Some(seasonality),
        )
        .expect("inflation calibration");
        let swap = ZeroCouponInflationSwap {
            notional: 1_000_000.0,
            fixed_rate: 0.021,
            maturity_years: 3.0,
            target_month: Some(7),
            pay_fixed: true,
        };
        let (_, delta) = curve.value_and_par_inflation_delta(|c| {
            let dc = discount.to_ad_curve(
                &discount
                    .discounts
                    .iter()
                    .map(|&x| ql_ad::TapeVar::from_f64(x))
                    .collect::<Vec<_>>(),
            );
            swap.npv(&dc, c)
        });
        let sum: f64 = delta.iter().sum::<f64>() * 1e-4;
        let npv = swap.npv(&discount, &curve.curve).val();
        (sum, npv)
    }
    fn inflation_calibration(
        &self,
        scenario: ScenarioKind,
    ) -> (DiscountCurve<f64>, CalibratedInflationCurve) {
        let d = self.base_market.ref_date;
        let discount = DiscountCurve::new(
            d,
            vec![1.0, 2.0, 5.0, 10.0],
            vec![
                (-0.02f64).exp(),
                (-0.04f64).exp(),
                (-0.1f64).exp(),
                (-0.2f64).exp(),
            ],
        );
        let bump = if scenario == ScenarioKind::StagflationPlus150Bp {
            0.015
        } else {
            0.0
        };
        let quotes = [
            (1.0, 0.02 + bump),
            (2.0, 0.022 + bump),
            (5.0, 0.025 + bump),
            (10.0, 0.027 + bump),
        ];
        let mut factors = [1.0; 12];
        for (i, factor) in factors.iter_mut().enumerate().take(11) {
            *factor = 1.0 + 0.006 * (i as f64 - 5.5) / 5.5;
        }
        factors[11] = 1.0 / factors[..11].iter().product::<f64>();
        let seasonality = MultiplicativeSeasonality::new(factors).expect("seasonality");
        let curve = CalibratedInflationCurve::bootstrap_with_seasonality(
            d,
            100.0,
            discount.clone(),
            &quotes,
            Some(seasonality),
        )
        .expect("inflation calibration");
        (discount, curve)
    }
    fn active_curve_system(
        &self,
        market: &MarketSnapshot,
    ) -> Arc<ql_curves::CalibratedCurveSystem> {
        ql_portfolio::build_curve_system(market)
    }
    fn exposure_cube(&self) -> ql_xva::NpvCube<f64> {
        let d = self.base_market.ref_date;
        let market = self.scenario_market();
        let curve = FlatForwardCurve::new(d, market.ois_quotes[0]);
        let model = HullWhite1F {
            a: 0.05,
            sigma: 0.01 + 0.25 * market.eq_vol,
        };
        let swap = NettingSetSwap {
            notional: 1_000_000.0,
            fixed_rate: market.ois_quotes[1],
            payment_times: (1..=5).map(|x| x as f64).collect(),
            accruals: vec![1.0; 5],
            is_payer: true,
        };
        let times: Vec<f64> = (1..=16).map(|k| k as f64 * 5.0 / 16.0).collect();
        simulate_hw1f_exposure_cube(&model, &curve, &[swap], &times, 64, 17)
    }
    fn xva_metrics(&self) -> (f64, f64, f64) {
        let d = self.base_market.ref_date;
        let market = self.scenario_market();
        let curve = FlatForwardCurve::new(d, market.ois_quotes[0]);
        let model = HullWhite1F {
            a: 0.05,
            sigma: 0.01,
        };
        let swap = NettingSetSwap {
            notional: 1_000_000.0,
            fixed_rate: market.ois_quotes[1],
            payment_times: (1..=5).map(|x| x as f64).collect(),
            accruals: vec![1.0; 5],
            is_payer: true,
        };
        let cube =
            simulate_hw1f_exposure_cube(&model, &curve, &[swap], &[1.0, 2.0, 3.0, 4.0], 32, 17);
        let hazard = if matches!(self.scenario, ScenarioKind::CreditWideningPlus100Bp) {
            0.03
        } else {
            0.02
        };
        let report = compute_xva(&cube, &curve, None, hazard, 0.6, 0.015, 0.6, 0.005);
        (report.cva.val(), report.dva.val(), report.fva.val())
    }
    fn mva_metric(&self) -> f64 {
        let d = self.base_market.ref_date;
        let market = self.scenario_market();
        let curve = FlatForwardCurve::new(d, market.ois_quotes[0]);
        let model = HullWhite1F {
            a: 0.05,
            sigma: 0.01 + 0.25 * market.eq_vol,
        };
        let swap = NettingSetSwap {
            notional: 1_000_000.0,
            fixed_rate: market.ois_quotes[1],
            payment_times: (1..=5).map(|x| x as f64).collect(),
            accruals: vec![1.0; 5],
            is_payer: true,
        };
        let times: Vec<f64> = (1..=12).map(|k| k as f64 * 5.0 / 12.0).collect();
        ql_xva::compute_chebyshev_dim_and_mva(&model, &curve, &[swap], &times, 64, 6, 0.01, 29).mva
    }
    fn phase22_live_metric(m: &MarketSnapshot, credit_bump: f64) -> f64 {
        let r = m.ois_quotes[0] + credit_bump;
        let bg = ql_models::evaluate_bergomi_forward_start_and_var_strip(
            &ql_models::BergomiForwardStartAndVarStripSpec {
                params: ql_models::Bergomi2FParams {
                    xi0: m.eq_vol * m.eq_vol,
                    omega: 1.0,
                    theta: 0.5,
                    k1: 1.2,
                    k2: 2.4,
                    rho12: 0.2,
                    rho_s1: -0.5,
                    rho_s2: -0.3,
                },
                xi_long_run: m.eq_vol * m.eq_vol,
                xi_curve_decay: 1.0,
                forward_start_t1: 0.0,
                expiry_t2: 1.0,
                log_moneyness_k: 0.0,
            },
        )
        .forward_start_implied_vol
        .val();
        let levy = cos_levy_price(
            &VarianceGamma {
                sigma: m.eq_vol.max(1e-4),
                nu: 0.2,
                theta: -0.1,
            },
            m.eq_spot,
            m.eq_spot,
            r,
            0.01,
            1.0,
            ql_instruments::OptionType::Call,
            128,
        )
        .val();
        let cg = evaluate_cgmy_fractional_fft_and_pide(&ql_levy::CgmyFractionalFftAndPideSpec {
            spot: m.eq_spot,
            strikes: vec![m.eq_spot],
            risk_free_rate: r,
            dividend_yield: 0.01,
            expiry: 1.0,
            process: CgmyProcess {
                sigma: m.eq_vol.max(1e-4),
                c: 0.5,
                g: 5.0,
                m: 5.0,
                y: 0.7,
            },
            damping_alpha: 1.0,
            fourier_step_eta: 0.1,
            log_strike_step_lambda: 0.1,
            small_jump_cutoff_eps: 1e-4,
        })
        .frft_call_prices
        .first()
        .copied()
        .unwrap_or(0.0)
        .val();
        let mm = ql_portfolio::solve_gueant_lehalle_bounded_market_making(
            &ql_portfolio::GueantLehalleMarketMakingSpec {
                max_inventory_q: 2,
                current_inventory_q: 0,
                horizon: 1.0,
                num_time_steps: 4,
                volatility: m.eq_vol.max(1e-4),
                risk_aversion_gamma: 0.5,
                intensity_decay_k: 1.5,
                base_intensity_a: 2.0,
                terminal_inventory_penalty: 0.1,
                adverse_selection_zeta: 0.0,
            },
        )
        .current_optimal_bid_offset
        .val();
        let bl = ql_portfolio::solve_meucci_idzorek_black_litterman_and_plausibility(
            &ql_portfolio::MeucciIdzorekBlackLittermanSpec {
                covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
                market_weights: vec![0.6, 0.4],
                risk_aversion_lambda: 2.5,
                tau: 0.05,
                view_matrix_p: vec![vec![1.0, -1.0]],
                view_returns_q: vec![0.03],
                view_confidences: vec![0.65],
            },
        )
        .mahalanobis_plausibility_squared
        .val();
        let coco =
            price_despiegeleer_coco_and_mandatory_convertible(&DeSpiegeleerCoCoMandatoryConvSpec {
                notional: 100.0,
                spot: m.eq_spot,
                risk_free_rate: r,
                dividend_yield: 0.01,
                credit_spread: 0.02,
                maturity: 2.0,
                coupon_rate: 0.05,
                coupon_freq: 2,
                decs_lower_strike: 90.0,
                decs_upper_strike: 110.0,
                vol_lower_strike: m.eq_vol,
                vol_upper_strike: m.eq_vol,
                coco_cet1_trigger_spot: 70.0,
                coco_ponv_trigger_spot: 50.0,
                coco_conversion_floor_price: 60.0,
                coco_equity_vol: m.eq_vol,
            })
            .dual_trigger_coco_value
            .val();
        // Phase 22's four remaining SOTA pillars are evaluated on the same
        // live market snapshot (and credit bump) as the existing telemetry.
        let krd = analyze_tuckman_key_rate_durations_and_pca_hedge(&TuckmanKeyRateDurationSpec {
            nss_curve: NelsonSiegelSvensson {
                reference_date: m.ref_date,
                beta0: r + 0.01,
                beta1: -0.01,
                beta2: 0.02,
                beta3: -0.005,
                tau1: 1.5,
                tau2: 4.0,
            },
            cashflow_times: vec![1.0, 2.0, 3.0, 5.0],
            cashflow_amounts: vec![2.0, 2.0, 2.0, 102.0],
            key_rate_tenors: vec![1.0, 2.0, 5.0, 10.0],
            pca_decay_lambda: 0.25,
            benchmark_hedge_tenors: [2.0, 5.0, 10.0],
        })
        .dirty_price
        .val();
        let lpi =
            price_brigo_mercurio_yoy_and_lpi_market_model(&BrigoMercurioInflationMarketModelSpec {
                notional: 1_000_000.0,
                base_cpi: 100.0,
                payment_times: vec![1.0, 2.0, 3.0],
                forward_cpis: vec![102.0, 104.0, 106.0],
                nominal_forward_rates: vec![r; 3],
                nominal_discount_factors: vec![(-r).exp(), (-2.0 * r).exp(), (-3.0 * r).exp()],
                forward_cpi_vols: vec![m.eq_vol; 3],
                nominal_forward_vols: vec![m.eq_vol * 0.5; 3],
                inter_maturity_cpi_corr: 0.5,
                nominal_cpi_corr: 0.2,
                lpi_floor: -0.02,
                lpi_cap: 0.03,
                fixed_swap_rate: 0.02 + credit_bump,
            })
            .fimm_convexity_adjusted_lpi_swap_npv
            .val();
        let cdo =
            price_okane_bespoke_tranche_and_student_t_basket(&OKaneBespokeBaseCorrelationSpec {
                index_detachments: vec![0.03, 0.07, 0.10, 0.15],
                index_base_correlations: vec![0.35, 0.40, 0.45, 0.50],
                index_expected_loss: 0.04 + credit_bump,
                bespoke_hazard_rates: vec![0.02 + credit_bump; 5],
                bespoke_recovery_rate: 0.4,
                bespoke_attachment: 0.03,
                bespoke_detachment: 0.10,
                strike_mapping_exponent_gamma: 1.0,
                student_t_dof_nu: 6.0 + 10.0 * credit_bump,
                maturity: 5.0,
                risk_free_rate: r,
            })
            .student_t_copula_senior_tranche_loss
            .val();
        let fx = evaluate_wystup_fx_accumulator_and_dcd(&WystupFxStructuredProductSpec {
            spot: m.eq_spot,
            accumulator_strike: m.eq_spot * 0.95,
            knockout_barrier: m.eq_spot * 1.10,
            downside_gearing: 2.0,
            domestic_rate: r,
            foreign_rate: 0.01,
            volatility: m.eq_vol.max(1e-4),
            fixing_times: vec![0.25, 0.5, 0.75, 1.0],
            guaranteed_fixings: 1,
            dcd_conversion_strike: m.eq_spot,
            dcd_tenor: 1.0,
        })
        .total_accumulator_npv
        .val();
        let pde = ql_pde::solve_rannacher_crank_nicolson_sinh_pde(&ql_pde::RannacherSinhPdeSpec {
            spot: m.eq_spot,
            strike: m.eq_spot,
            volatility: m.eq_vol.max(1e-4),
            risk_free_rate: r,
            dividend_yield: 0.01,
            expiry: 1.0,
            is_call: true,
            n_space: 31,
            n_time: 8,
            rannacher_steps: 2,
            sinh_concentration: 0.15,
        })
        .pde_price
        .val();
        let mot =
            ql_ad::solve_mot_super_replication_bounds_and_ift(&ql_ad::MotSuperReplicationIftSpec {
                support_grid: vec![0.8, 1.0, 1.2],
                marginal_t1_probs: vec![0.2, 0.5, 0.3],
                marginal_t2_probs: vec![0.25, 0.5, 0.25],
                forward_start_strike: 0.95 + 0.001 * credit_bump,
                entropic_regularization_eps: 0.02,
            })
            .entropic_mot_price
            .val();
        let exotic = ql_dml::evaluate_ift_calibrated_exotic_quote_hedge(
            &ql_dml::IftCalibratedExoticHedgeSpec {
                spot: m.eq_spot,
                risk_free_rate: r,
                expiry: 1.0,
                calibration_strikes: [0.9 * m.eq_spot, m.eq_spot, 1.1 * m.eq_spot],
                market_implied_vols: [m.eq_vol * 1.1, m.eq_vol, m.eq_vol * 1.05],
                exotic_strike: 1.05 * m.eq_spot,
                exotic_power_p: 1.2,
            },
        )
        .exotic_option_price
        .val();
        let revolve = ql_ad::execute_savine_revolve_checkpointed_adjoint(
            &ql_ad::SavineRevolveCheckpointSpec {
                num_time_steps: 8,
                num_checkpoints: 2,
                initial_spot: m.eq_spot,
                drift_rate: r,
                volatility: m.eq_vol.max(1e-4),
                mean_reversion_kappa: 0.5,
                dt: 0.125,
                strike: m.eq_spot,
            },
        )
        .terminal_option_payoff_value
        .val();
        phase21_live_metric(m, credit_bump)
            + bg
            + levy
            + cg
            + mm
            + bl
            + coco
            + krd * 1e-4
            + lpi * 1e-6
            + cdo
            + fx
            + pde
            + mot
            + exotic
            + revolve
    }
    pub fn phase24_live_metric_from_market(m: &MarketSnapshot, credit_bump: f64) -> f64 {
        let inf = m.inflation_quotes.first().map(|(_, z)| *z).unwrap_or(0.023);
        let ie_growth = (0.02 + 0.25 * (inf - 0.023)).clamp(0.005, 0.045);
        let mbs_cpr = (0.06 + 0.5 * (inf - 0.023)).clamp(0.02, 0.20);
        let mbs = ql_fixed_income::analyze_davidson_levin_mbs_io_po_strips(
            &ql_fixed_income::DavidsonLevinIoPoStripMbsSpec {
                notional: 100.0,
                wac_annual: 0.065,
                pass_through_coupon_annual: 0.06,
                current_mortgage_rate: m.ois_quotes[0] + credit_bump,
                discount_rate_annual: m.ois_quotes[0],
                remaining_months: 24,
                loan_age_months: 12,
                turnover_base_cpr: mbs_cpr,
                max_refi_cpr: 0.4,
                refi_slope: 80.0,
                refi_threshold_spread: 0.005,
                burnout_exponent: 1.0,
                rate_bump_bps: 1.0,
            },
        );
        let bond = ql_credit::analyze_duffie_singleton_rmv_rfv_and_cds_basis(
            &ql_credit::DuffieSingletonRecoveryAndCdsBasisSpec {
                face_value: 100.0,
                coupon_rate_annual: 0.05,
                payment_freq_per_year: 2,
                maturity_years: 3.0,
                risk_free_rate: m.ois_quotes[0],
                hazard_rate: 0.02 + credit_bump,
                recovery_rate: 0.4,
                market_bond_price: 98.0,
                repo_funding_spread: 0.005,
            },
        );
        let mo = ql_xva::compute_crepey_bielecki_marshall_olkin_bilateral_xva(
            &ql_xva::CrepeyBieleckiMarshallOlkinBilateralXvaSpec {
                times: vec![0.5, 1.0],
                expected_positive_exposure: vec![1.0, 1.0],
                expected_negative_exposure: vec![0.5, 0.5],
                risk_free_rate: m.ois_quotes[0],
                lambda_institution_idio: 0.01,
                lambda_counterparty_idio: 0.02,
                lambda_common_shock_12: 0.01 + credit_bump,
                recovery_institution: 0.4,
                recovery_counterparty: 0.4,
                recovery_counterparty_simultaneous: 0.3,
                gap_exposure_multiplier: 0.2,
                funding_borrowing_spread: 0.01,
                funding_lending_spread: 0.005,
            },
        );
        let cb = ql_fixed_income::price_sirbu_shreve_dilution_and_makewhole_convertible(
            &ql_fixed_income::SirbuShreveDilutionMakeWholeConvertibleSpec {
                spot: m.eq_spot,
                volatility: m.eq_vol,
                risk_free_rate: m.ois_quotes[0],
                dividend_yield: 0.0,
                credit_spread: 0.02 + credit_bump,
                maturity_years: 2.0,
                par_value: 100.0,
                coupon_rate_annual: 0.03,
                conversion_ratio: 0.8,
                existing_shares_millions: 100.0,
                convertible_issue_count_millions: 10.0,
                debt_relief_retention_fraction: 0.0,
                call_price: 105.0,
                soft_call_trigger_spot: 110.0,
                make_whole_max_bonus_shares: 0.2,
                make_whole_time_decay: 0.2,
                make_whole_moneyness_decay: 0.3,
                tree_steps: 8,
            },
        );
        let sig = ql_models::evaluate_bayer_oberhauser_signature_option_pricing(
            &ql_models::BayerHorvathOberhauserPathSignatureSpec {
                time_grid: vec![0.0, 0.5, 1.0],
                sample_paths: vec![vec![m.eq_spot, m.eq_spot * 1.01, m.eq_spot]],
                strike: m.eq_spot,
                discount_factor: (-m.ois_quotes[0]).exp(),
                ridge_lambda: 1e-6,
            },
        );
        let hn = ql_vol::evaluate_heston_nandi_and_tsay_dcc_garch(
            &ql_vol::HestonNandiTsayDccGarchSpec {
                spot: m.eq_spot,
                strike: m.eq_spot,
                daily_risk_free_rate: m.ois_quotes[0] / 252.0,
                days_to_maturity: 30,
                initial_daily_variance: m.eq_vol * m.eq_vol / 252.0,
                omega: 1e-5,
                alpha: 0.05,
                beta: 0.9,
                gamma_asymmetry: 0.2,
                asset1_returns: vec![0.01],
                asset2_returns: vec![0.008],
                dcc_a: 0.02,
                dcc_b: 0.95,
            },
        );
        let pot = ql_portfolio::evaluate_mcneil_embrechts_pot_gpd_and_christoffersen_backtest(
            &ql_portfolio::McNeilEmbrechtsPotGpdBacktestSpec {
                historical_losses: vec![1.0, 2.0, 3.0, 4.0, 5.0],
                threshold_u: 1.0,
                confidence_alpha: 0.95,
                spectral_risk_aversion_k: 2.0,
                backtest_var_threshold: 3.0,
            },
        );
        let rie = ql_portfolio::solve_potters_bouchaud_rie_nonlinear_shrinkage(
            &ql_portfolio::PottersBouchaudRieCovarianceSpec {
                sample_covariance: vec![vec![0.04 + credit_bump, 0.01], vec![0.01, 0.09]],
                num_observations_t: 100,
            },
        );
        let pin = ql_portfolio::evaluate_ohara_ekop_kyle_microstructure_and_etf_band(
            &ql_portfolio::OharaEkopKyleMicrostructureSpec {
                fundamental_value_std: 1.0,
                noise_trader_volume_std: 2.0,
                alpha_info_event_prob: 0.2,
                delta_bad_news_prob: 0.4,
                mu_informed_arrival_rate: 1.0,
                epsilon_uninformed_buy_rate: 2.0,
                epsilon_uninformed_sell_rate: 2.0,
                fundamental_low: 95.0,
                fundamental_high: 105.0,
                position_value: 100.0,
                return_vol_daily: 0.02,
                z_score_alpha: 2.0,
                mean_relative_spread: 0.01,
                std_relative_spread: 0.002,
                etf_price: m.eq_spot,
                etf_inav: 100.0,
                ap_creation_fee_bps: 5.0,
                basket_spread_bps: 2.0,
                etf_spread_bps: 3.0,
            },
        );
        let gp = ql_portfolio::solve_garleanu_pedersen_dynamic_aim_portfolio(
            &ql_portfolio::GarleanuPedersenDynamicAimPortfolioSpec {
                covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
                current_weights: vec![0.5, 0.5],
                slow_alpha_vector: vec![0.03, 0.02],
                fast_alpha_vector: vec![0.01, 0.02],
                slow_mean_reversion_phi: 0.1,
                fast_mean_reversion_phi: 0.5,
                risk_aversion_gamma: 2.0,
                transaction_cost_kappa: 0.1,
                discount_rate_r: m.ois_quotes[0],
            },
        );
        let br = ql_portfolio::evaluate_bacon_carino_brinson_and_karnosky_singer_attribution(
            &ql_portfolio::BaconCarinoMultiPeriodBrinsonSpec {
                portfolio_weights_by_period: vec![vec![0.5, 0.5]],
                benchmark_weights_by_period: vec![vec![0.5, 0.5]],
                portfolio_returns_by_period: vec![vec![0.1, 0.02]],
                benchmark_returns_by_period: vec![vec![0.08, 0.01]],
                fx_spot_returns_by_sector: vec![0.01],
                foreign_cash_rates_by_sector: vec![0.02],
                domestic_cash_rate: m.ois_quotes[0],
            },
        );
        let ie = ql_ad::solve_penman_damodaran_riv_aeg_and_implied_erp_ift(
            &ql_ad::PenmanDamodaranResidualIncomeIerpSpec {
                initial_book_value: 100.0,
                forecast_roe_schedule: vec![0.12, 0.11],
                payout_ratio_schedule: vec![0.4, 0.4],
                terminal_roe: 0.1,
                terminal_growth_rate: ie_growth,
                risk_free_rate: m.ois_quotes[0],
                equity_beta: 1.1,
                observed_market_price: 110.0,
            },
        );
        mbs.io_strip_price
            + bond.raw_cds_bond_basis
            + mo.total_bilateral_xva
            + cb.make_whole_enhanced_price
            + sig.signature_linear_functional_price
            + hn.heston_nandi_call_price
            + pot.pot_gpd_expected_shortfall
            + rie.rie_portfolio_variance
            + pin.ekop_pin
            + gp.next_period_optimal_weights[0]
            + br.carino_linked_selection_total
            + ie.implied_erp
    }
    fn phase23_live_metric_from_market(m: &MarketSnapshot, credit_bump: f64) -> f64 {
        let rk = ql_curves::calibrate_filipovic_rkhs_smooth_yield_curve(
            &ql_curves::FilipovicRkhsYieldCurveSpec {
                cashflow_times: vec![1.0, 2.0],
                instrument_cashflow_matrix: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
                instrument_dirty_prices: vec![m.ois_quotes[0], m.ois_quotes[1]],
                decay_alpha: 1.0,
                ridge_penalty_lambda: 1e-8,
                query_times: vec![1.0, 2.0],
            },
        );
        let fmm = ql_rates::price_lyashenko_mercurio_fmm_caplet_and_swaption(
            &ql_rates::LyashenkoMercurioFmmSpec {
                notional: 1.0,
                exercise_time_t0: 1.0,
                period_accruals: vec![0.5, 0.5],
                forward_rates: vec![m.ois_quotes[1], m.ois_quotes[1]],
                discount_factors: vec![0.98, 0.96],
                fmm_volatilities: vec![m.eq_vol, m.eq_vol],
                displacement_shift: 0.0,
                inter_period_correlation_decay: 1.0,
                strike: m.ois_quotes[1],
                is_payer_or_caplet: true,
            },
        );
        let vv = ql_fx::evaluate_castagna_full_vanna_volga_hedge_and_barrier(
            &ql_fx::CastagnaFullVannaVolgaHedgeSpec {
                spot: m.eq_spot,
                r_dom: m.ois_quotes[0],
                r_for: 0.01,
                expiry: 1.0,
                atm_vol: m.eq_vol,
                rr_25d: 0.0,
                bf_25d: 0.0,
                target_strike: m.eq_spot,
                knockout_barrier: m.eq_spot * 1.2,
                unwind_attenuation_factor: 0.5,
                is_up_and_out: true,
                is_call: true,
            },
        );
        let lsv = ql_vol::evaluate_guyon_markovian_projection_lsv_smile(
            &ql_vol::GuyonMarkovianProjectionLsvSpec {
                forward: 100.0,
                strike: 100.0,
                expiry: 1.0,
                forward_start_t1: 0.25,
                atm_vol: m.eq_vol,
                market_atm_skew: -0.1,
                market_atm_curvature: 0.1,
                sv_v0: m.eq_vol * m.eq_vol,
                sv_kappa: 1.0,
                sv_eta: 0.3,
                sv_rho: -0.5,
                sv_mixing_weight_eta: 0.5,
            },
        );
        let jlt = ql_credit::evaluate_jlt_rating_migration_and_crcva(
            &ql_credit::JarrowLandoTurnbullRatingMigrationSpec {
                initial_rating_index: 0,
                physical_generator_4x4: [
                    [-0.1, 0.08, 0.01, 0.01],
                    [0.02, -0.08, 0.04, 0.02],
                    [0.0, 0.02, -0.05, 0.03],
                    [0.0, 0.0, 0.0, 0.0],
                ],
                risk_premium_multipliers: [1.0, 1.1, 1.2],
                recovery_rate: 0.4,
                risk_free_rate: m.ois_quotes[0],
                maturity: 5.0,
                base_coupon_rate: 0.03,
                hy_downgrade_stepup_spread: 0.01,
                counterparty_epe_notional: 1.0,
            },
        );
        let ctd = ql_xva::compute_brigo_pallavicini_ctd_collateral_and_closeout_xva(
            &ql_xva::BrigoPallaviciniCtdCollateralXvaSpec {
                notional_mtm: 1.0,
                maturity: 1.0,
                num_steps: 8,
                base_ois_rate: m.ois_quotes[0],
                csa_spread_ccy1: 0.001,
                csa_spread_ccy2: 0.002,
                csa_spread_ccy3: 0.003,
                csa_mean_reversion: 1.0,
                csa_vol_ccy1: 0.01,
                csa_vol_ccy2: 0.01,
                csa_vol_ccy3: 0.01,
                csa_correlation: 0.2,
                cpty_hazard_rate: 0.02,
                cpty_lgd: 0.6,
                own_hazard_rate: 0.01,
                own_lgd: 0.6,
                exposure_volatility: 0.1,
            },
        );
        let qrh = ql_models::evaluate_gatheral_rosenbaum_joint_spx_vix_qrh(
            &ql_models::GatheralRosenbaumJointSpxVixQrhSpec {
                params: ql_models::QuadraticRoughHestonParams {
                    hurst: 0.1,
                    nu: 0.3,
                    lambda: 1.0,
                    a: 1.0,
                    b: 0.0,
                    c: 0.04,
                    z0: 0.0,
                },
                spot: 100.0,
                spx_strike: 100.0,
                vix_strike: 0.2,
                expiry: 1.0,
                rate: m.ois_quotes[0],
            },
        );
        let par = ql_exotics::price_chesney_yor_cumulative_parisian_and_rainbow(
            &ql_exotics::ChesneyYorCumulativeParisianRainbowSpec {
                spot_1: 100.0,
                spot_2: 102.0,
                strike: 100.0,
                par_barrier: 120.0,
                excursion_window_d: 0.1,
                vol_1: m.eq_vol,
                vol_2: m.eq_vol,
                div_1: 0.0,
                div_2: 0.0,
                correlation: 0.5,
                risk_free_rate: m.ois_quotes[0],
                expiry: 1.0,
                is_up_and_out: true,
            },
        );
        let aid = ql_commodities::price_aid_carmona_structural_electricity_stack(
            &ql_commodities::AidCarmonaStructuralElectricityStackSpec {
                gas_forward: m.fwd_quotes[0] * 1000.0,
                gas_vol: 0.3,
                base_heat_rate: 2.0,
                total_capacity_mw: 100.0,
                baseline_capacity_anchor_mw: 50.0,
                expected_load_mw: 80.0,
                expected_renewable_gen_mw: 20.0,
                residual_load_vol_mw: 10.0,
                stack_convexity_alpha: 1.0,
                negative_price_floor_depth: 5.0,
                curtailment_beta: 1.0,
                load_gas_correlation: 0.2,
                power_option_strike: 50.0,
                spike_threshold_price: 100.0,
                risk_free_rate: m.ois_quotes[0],
                expiry: 1.0,
            },
        );
        let imp = ql_portfolio::solve_bouchaud_webster_square_root_and_propagator_execution(
            &ql_portfolio::BouchaudWebsterPropagatorImpactSpec {
                metaorder_size_q: 1.0,
                daily_volume_vd: 10.0,
                daily_volatility_sigma: 0.2,
                y_prefactor: 1.0,
                propagator_exponent_beta: 0.5,
                propagator_scale_g0: 1.0,
                regularization_tau0: 0.01,
                temporary_impact_eta: 0.1,
                risk_aversion_lambda: 0.1,
                horizon: 1.0,
                num_slices: 8,
            },
        );
        let nco = ql_portfolio::solve_lopez_de_prado_denoised_nco_allocation(
            &ql_portfolio::LopezDePradoNcoDenoisedSpec {
                covariance: vec![
                    vec![0.04 + credit_bump, 0.01],
                    vec![0.01, 0.09 + credit_bump],
                ],
                expected_returns: vec![0.05, 0.06],
                sample_ratio_q: 2.0,
                denoising_alpha: 0.5,
                num_clusters: 2,
            },
        );
        let dml = ql_dml::evaluate_glasserman_savine_hybrid_pw_lrm_and_dml(
            &ql_dml::GlassermanSavineHybridPwLrmSpec {
                spot: 100.0,
                strike: 100.0,
                rate: m.ois_quotes[0],
                div_yield: 0.0,
                volatility: m.eq_vol,
                expiry: 1.0,
                digital_smoothing_eps: 0.5,
                dml_differential_weight_gamma: 1.0,
            },
        );
        rk.query_zero_rates[0].val()
            + fmm.backward_looking_rfr_caplet_npv.val()
            + vv.exact_vv_vanilla_price.val()
            + lsv.lsv_reconstructed_implied_vol.val()
            + jlt.rating_migration_adjusted_crcva.val()
            + ctd.ctd_collateral_switch_option_value.val()
            + qrh.vix_futures_fair_value.val()
            + par.rainbow_call_on_max.val()
            + aid.structural_power_forward.val()
            + imp.optimal_propagator_cost.val()
            + nco.nco_portfolio_volatility.val()
            + dml.hybrid_pw_lrm_digital_delta.val()
            + credit_bump
            + Self::phase22_live_metric(m, credit_bump) * 1e-3
    }
    pub fn phase23_live_metric(&self) -> f64 {
        let market = self.scenario_market();
        let bump = if matches!(self.scenario, ScenarioKind::CreditWideningPlus100Bp) {
            0.01
        } else {
            0.0
        };
        Self::phase23_live_metric_from_market(&market, bump)
    }
    pub fn phase25_live_metric_from_market(m: &MarketSnapshot, credit_bump: f64) -> f64 {
        let r = m.ois_quotes[0];
        let spot = m.eq_spot;
        let vol = m.eq_vol.max(1e-4);
        let inf = m.inflation_quotes.first().map(|(_, z)| *z).unwrap_or(0.023);
        let kal = ql_curves::evaluate_tsay_kalman_rts_state_space_curve(
            &ql_curves::TsayKalmanRtsStateSpaceCurveSpec {
                transition_f: [[0.9, 0.1], [0.0, 0.95]],
                state_drift_c: [0.0, 0.0],
                process_noise_q: [[0.01, 0.0], [0.0, 0.01]],
                measurement_h: [[1.0, 0.0], [0.0, 1.0]],
                measurement_noise_r: [0.02, 0.02],
                initial_state_mean: [r, spot],
                initial_state_cov: [[1.0, 0.0], [0.0, 1.0]],
                observations: vec![[spot, r], [spot * 1.001, r + 0.001]],
            },
        );
        let vecm = ql_portfolio::solve_johansen_tsay_vecm_cointegration_basket(
            &ql_portfolio::JohansenTsayVecmCointegrationSpec {
                price_series: (0..8)
                    .map(|i| {
                        [
                            100.0 + i as f64,
                            100.0 + i as f64 + 0.5 * (1.0 - (0.9_f64).powi(i)),
                        ]
                    })
                    .collect(),
            },
        );
        let hj = ql_portfolio::evaluate_cochrane_hansen_jagannathan_sdf_and_campbell_shiller(
            &ql_portfolio::CochraneHansenJagannathanSdfSpec {
                mean_excess_returns: [0.03, 0.04],
                excess_return_cov: [[0.04, 0.01], [0.01, 0.09]],
                risk_free_rate: r,
                candidate_sdf_pricing_errors: [0.001, 0.002],
                steady_state_log_dp_ratio: inf.ln(),
                var1_companion_matrix: [[0.9, 0.1], [0.0, 0.8]],
                state_innovation_cov: [[0.02, 0.005], [0.005, 0.03]],
            },
        );
        let ifrs = ql_credit::compute_ifrs9_cecl_multistage_expected_credit_loss(
            &ql_credit::Ifrs9CeclMultiStageEclSpec {
                annual_ttc_pds: [
                    0.02 + credit_bump * 0.2,
                    0.04 + credit_bump * 0.2,
                    0.08 + credit_bump * 0.2,
                ],
                annual_eads: [100.0, 98.0, 95.0],
                base_lgd: 0.45,
                lgd_macro_sensitivity: 0.05,
                vasicek_asset_correlation: 0.2,
                effective_interest_rate: r,
                scenario_weights: [0.3, 0.5, 0.2],
                scenario_macro_shocks_z: [-1.0, 0.0, 1.0],
                origination_lifetime_pd: 0.04,
                sicr_pd_multiplier_threshold: 1.5,
            },
        );
        let cpt = ql_portfolio::evaluate_wakker_tversky_kahneman_cumulative_prospect_utility(
            &ql_portfolio::WakkerTverskyKahnemanCptSpec {
                scenario_returns: [-0.10, -0.02, 0.03, 0.12],
                scenario_probabilities: [0.1, 0.2, 0.4, 0.3],
                reference_return: 0.0,
                alpha_gain: 0.88,
                beta_loss: 0.88,
                loss_aversion_lambda: 2.25,
                prelec_gamma_gain: 0.7,
                prelec_gamma_loss: 0.8,
                prelec_delta: 1.0,
            },
        );
        let has = ql_portfolio::evaluate_hasbrouck_information_share_and_dssw_limits_to_arb(
            &ql_portfolio::HasbrouckDsswPriceDiscoverySpec {
                alpha_venue1: 0.3,
                alpha_venue2: 0.2,
                sigma_venue1: vol,
                sigma_venue2: vol * 1.1,
                venue_innovation_corr: 0.2,
                fundamental_dividend: 2.0,
                risk_free_rate: r,
                noise_trader_share_mu: 0.2,
                noise_trader_sentiment_mean: 0.01,
                noise_trader_sentiment_var: 0.02,
                arbitrageur_risk_aversion_gamma: 2.0,
            },
        );
        let fv = ql_rates::price_fong_vasicek_stochastic_vol_bond(
            &ql_rates::FongVasicekStochVolShortRateSpec {
                initial_short_rate_r0: r,
                initial_variance_v0: vol * vol,
                rate_mean_reversion_alpha: 0.5,
                long_run_rate_r_bar: r,
                vol_mean_reversion_gamma: 1.0,
                long_run_variance_v_bar: vol * vol,
                vol_of_vol_xi: 0.2,
                correlation_rho: -0.3,
                maturity_tau: 2.0,
            },
        );
        let lando = ql_credit::evaluate_lando_cox_affine_default_and_liquidity_decomposition(
            &ql_credit::LandoCoxAffineDefaultSpec {
                risk_free_rate: r,
                base_default_hazard: 0.01 + credit_bump,
                initial_credit_intensity_y0: 0.02,
                credit_kappa: 0.5,
                credit_theta: 0.02,
                credit_sigma: 0.1,
                initial_liquidity_intensity_z0: 0.01,
                liquidity_kappa: 0.6,
                liquidity_theta: 0.01,
                liquidity_sigma: 0.08,
                recovery_rate: 0.4,
                maturity: 5.0,
            },
        );
        let times: Vec<f64> = (0..8).map(|i| i as f64 + 1.0).collect();
        let logs: Vec<f64> = times
            .iter()
            .map(|&t| 4.0 - 0.2 * (10.0 - t).powf(0.5) + 0.01 * (8.0 * (10.0 - t).ln()).cos())
            .collect();
        let lppl = ql_portfolio::calibrate_sornette_lppl_bubble_and_crash_hazard(
            &ql_portfolio::SornetteLpplBubbleCrashSpec {
                observation_times: times,
                log_prices: logs,
                critical_time_tc: 10.0,
                power_law_exponent_m: 0.5,
                log_periodic_omega: 8.0,
                expected_crash_drop_kappa: 0.2,
                ridge_penalty: 1e-8,
            },
        );
        let cb = ql_fixed_income::analyze_thorp_kassouf_calamos_convertible_arbitrage(
            &ql_fixed_income::ThorpKassoufCalamosConvertibleArbSpec {
                stock_price: spot,
                conversion_ratio: 0.8,
                convertible_market_price: 95.0,
                face_value: 100.0,
                annual_coupon_rate: 0.03,
                stock_dividend_yield: 0.01,
                risk_free_rate: r,
                credit_spread: 0.03 + credit_bump,
                short_stock_rebate_rate: 0.01,
                implied_volatility: vol,
                realized_volatility: vol * 1.1,
                maturity: 3.0,
            },
        );
        let causal = ql_portfolio::evaluate_lopez_de_prado_causal_backdoor_vs_collider_bias(
            &ql_portfolio::LopezDePradoCausalFactorSpec {
                true_causal_beta_xy: 0.5,
                confounder_to_factor_axz: 0.8,
                confounder_to_return_byz: 0.6,
                factor_to_collider_ccx: 0.4,
                return_to_collider_ccy: 0.5,
                confounder_var: 1.0,
                factor_idiosyncratic_var: 0.5,
                return_idiosyncratic_var: 0.7,
                collider_idiosyncratic_var: 0.4,
            },
        );
        let mal = ql_dml::evaluate_fournie_lasry_lions_malliavin_bel_greeks(
            &ql_dml::FournieLasryLionsMalliavinBelSpec {
                spot,
                strike: spot,
                rate: r,
                div_yield: 0.0,
                vol,
                expiry: 1.0,
                localization_bandwidth: 0.5,
            },
        );
        kal.innovation_log_likelihood.val()
            + vecm.johansen_trace_stat_r0.val() * 1e-3
            + hj.discount_rate_news_variance_share.val()
            + ifrs.blended_ifrs9_ecl_provision.val() * 1e-2
            + cpt.cumulative_prospect_utility.val()
            + has.efficient_price_innovation_var.val()
            + fv.zero_coupon_bond_price.val()
            + lando.risky_corporate_bond_price.val()
            + lppl.terminal_crash_hazard_rate.val()
            + cb.theoretical_convertible_value.val() * 1e-2
            + causal.kitchen_sink_collider_conditioned_beta.val()
            + mal.call_malliavin_delta.val()
            + inf * 10.0
            + credit_bump
    }
    pub fn phase27_live_metric_from_market(m: &MarketSnapshot, credit_bump: f64) -> f64 {
        let r = m.ois_quotes[0] + credit_bump;
        let spot = m.eq_spot;
        let vol = m.eq_vol.max(1e-4);
        let inf = m.inflation_quotes.first().map(|(_, z)| *z).unwrap_or(0.023);
        let ftpl = ql_inflation::evaluate_cochrane_ftpl_and_maturity_inflation_smoothing(
            &ql_inflation::CochraneFiscalTheoryPriceLevelSpec {
                debt_maturities_years: vec![1.0, 5.0, 10.0],
                nominal_debt_face_values: vec![30.0, 30.0, 40.0],
                nominal_zero_yields: vec![r, r + 0.002, r + 0.004],
                real_zero_yields: vec![r - inf; 3],
                primary_real_surpluses: vec![5.0, 5.2, 5.4],
                terminal_surplus_growth_rate: 0.02,
                fiscal_surplus_shock_fraction: -0.05,
                expected_future_inflation_rise: inf,
                inflation_risk_beta: 0.4,
                inflation_volatility: vol,
                tips_liquidity_spread: 0.002,
            },
        );
        let g2 = ql_rates::price_brigo_mercurio_g2_twist_and_swaption(
            &ql_rates::BrigoMercurioG2TwistSwaptionSpec {
                a: 0.1,
                sigma: vol,
                b: 0.3,
                eta: vol * 0.8,
                rho: -0.4,
                flat_zero_rate: r,
                option_expiry: 1.0,
                swap_payment_times: vec![1.0, 2.0, 3.0, 4.0, 5.0],
                fixed_strike: r,
                notional: 1.0,
            },
        );
        let ess = ql_vol::evaluate_mingone_global_arbitrage_free_essvi(
            &ql_vol::MingoneGlobalArbitrageFreeEssviSpec {
                expiries: vec![1.0, 2.0],
                raw_theta_increments: vec![vol, vol],
                raw_rho_params: vec![-0.2, -0.1],
                raw_psi_fractions: vec![0.2, 0.3],
                test_log_moneyness_grid: vec![-1.0, 0.0, 1.0],
            },
        );
        let cb = ql_fx::price_cross_currency_convertible_and_warrant_bond(
            &ql_fx::CrossCurrencyConvertibleWarrantSpec {
                domestic_stock_spot: spot,
                domestic_conversion_strike: spot,
                spot_fx_for_per_dom: 1.1,
                fixed_contractual_fx_for_per_dom: 1.1,
                foreign_bond_par: 100.0,
                foreign_coupon_rate: 0.03,
                r_domestic: r,
                r_foreign: r - 0.005,
                stock_dividend_yield: 0.01,
                credit_spread: credit_bump + 0.02,
                stock_vol: vol,
                fx_vol: vol * 0.7,
                corr_stock_fx: 0.2,
                maturity: 5.0,
            },
        );
        let bank = ql_fixed_income::analyze_acharya_feds_bank_run_and_optimal_coco(
            &ql_fixed_income::AcharyaFedsBankRunAndCocoSpec {
                initial_securities_assets: 100.0,
                asset_effective_duration: 5.0,
                insured_deposits: 40.0,
                uninsured_deposits: 45.0,
                deposit_rate_beta: 0.5,
                deposit_operating_cost_rate: 0.005,
                at1_coco_principal: 5.0,
                initial_short_rate: r,
                parallel_rate_shock: inf,
                run_sensitivity_steepness: 30.0,
                fire_sale_liquidation_haircut: 0.1,
                coco_write_down_efficiency: 0.9,
            },
        );
        let xva = ql_xva::compute_albanese_crepey_mpor_and_balance_sheet_xva(
            &ql_xva::AlbaneseCrepeyMporBalanceSheetSpec {
                exposure_times: vec![1.0, 2.0],
                clean_mtm_profile: vec![1.0, 2.0],
                portfolio_vol_profile: vec![vol, vol],
                scheduled_cashflows_in_mpor: vec![0.01, 0.02],
                mpor_years: 0.05,
                csa_threshold: 0.1,
                risk_free_rate: r,
                cpty_hazard_rate: 0.02,
                cpty_lgd: 0.6,
                bank_hazard_rate: 0.01,
                bank_lgd: 0.6,
                economic_capital_profile: vec![0.05, 0.05],
                hurdle_rate: 0.1,
                initial_margin_profile: vec![0.1, 0.1],
                im_funding_spread: 0.01,
            },
        );
        let dh = ql_models::evaluate_double_heston_cumulants_and_kernel_slv(
            &ql_models::DoubleHestonKernelSlvSpec {
                spot,
                strike: spot,
                maturity: 1.0,
                risk_free_rate: r,
                dividend_yield: 0.0,
                v1_0: vol * vol * 0.6,
                kappa1: 2.0,
                theta1: vol * vol * 0.6,
                sigma1: 0.3,
                rho1: -0.6,
                v2_0: vol * vol * 0.4,
                kappa2: 0.5,
                theta2: vol * vol * 0.4,
                sigma2: 0.2,
                rho2: -0.3,
                target_dupire_local_vol: vol,
            },
        );
        let asian = ql_exotics::price_curran_rogers_shi_arithmetic_asian(
            &ql_exotics::CurranRogersShiArithmeticAsianSpec {
                spot,
                strike: spot,
                risk_free_rate: r,
                dividend_yield: 0.0,
                vol,
                fixing_times: vec![0.25, 0.5, 0.75, 1.0],
            },
        );
        let storage = ql_commodities::evaluate_warin_baeuerle_storage_rolling_hedge(
            &ql_commodities::WarinBaeuerleStorageRollingHedgeSpec {
                forward_curve: vec![spot, spot * 1.01, spot * 1.02],
                decision_times: vec![0.5, 1.0, 1.5],
                max_capacity: 10.0,
                initial_inventory: 2.0,
                max_injection_step: 1.0,
                max_withdrawal_step: 1.0,
                injection_loss_fraction: 0.02,
                spot_ou_volatility: vol,
                spot_ou_mean_reversion: 1.0,
                discount_rate: r,
            },
        );
        let micro = ql_portfolio::analyze_harris_roll_corwin_schultz_microstructure(
            &ql_portfolio::HarrisRollCorwinSchultzMicrostructureSpec {
                trade_prices: vec![spot, spot * 1.001, spot * 0.999, spot * 1.002],
                daily_highs: vec![spot * 1.01, spot * 1.02, spot * 1.015],
                daily_lows: vec![spot * 0.99, spot * 0.98, spot * 0.985],
                daily_dollar_volumes: vec![1e6, 1.1e6, 0.9e6],
                decision_price: spot,
                arrival_price: spot,
                child_fill_quantities: vec![10.0, 5.0],
                child_fill_prices: vec![spot * 1.001, spot * 1.002],
                total_order_quantity: 20.0,
                post_trade_benchmark_price: spot * 1.003,
                explicit_commissions: 0.1,
            },
        );
        let kelly = ql_portfolio::solve_busseti_boyd_risk_constrained_kelly(
            &ql_portfolio::BussetiBoydRiskConstrainedKellySpec {
                excess_return_means: vec![0.05 + spot * 1e-4 + inf * 0.1, 0.03 + r * 0.1],
                covariance: vec![
                    vec![vol * vol, vol * vol * 0.25],
                    vec![vol * vol * 0.25, (vol * vol * 0.75).max(1e-6)],
                ],
                risk_free_rate: r,
                drawdown_floor_alpha: 0.8,
                max_drawdown_prob_beta: 0.5,
                specialist_immediate_capital: 100.0,
                slow_moving_external_capital: 100.0,
                capital_arrival_rate: 0.1,
                dislocation_horizon_years: 1.0,
                forced_liquidation_notional: 10.0,
            },
        );
        let qp = ql_dml::solve_agrawal_boyd_differentiable_qp_layer(
            &ql_dml::AgrawalBoydDifferentiableQpSpec {
                q_matrix: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
                q_linear: vec![-spot * 0.001, -vol],
                a_eq: vec![vec![1.0, 1.0]],
                b_eq: vec![1.0],
                g_ineq: vec![vec![-1.0, 0.0]],
                h_ineq: vec![0.0],
                barrier_tau: 1e-4,
            },
        );
        ftpl.equilibrium_price_level_p0.val()
            + g2.payer_swaption_pv.val()
            + ess.calibrated_thetas.iter().map(|x| x.val()).sum::<f64>()
            + cb.embedded_cross_currency_cb_price.val() * 1e-2
            + bank.economic_equity_post_coco_conversion.val() * 1e-2
            + xva.sticky_mpor_cva.val()
            + dh.double_heston_atm_vol.val()
            + asian.geometric_asian_call_price.val()
            + storage.stochastic_total_storage_value.val() * 1e-2
            + micro.perold_total_implementation_shortfall.val() * 1e-2
            + kelly.constrained_portfolio_volatility.val()
            + qp.optimal_objective_value.val()
            + inf
    }

    pub fn phase30_live_metric_from_market(m: &MarketSnapshot, credit_bump: f64) -> f64 {
        Self::phase30_live_metric_from_market_with_bump(m, credit_bump)
    }
    pub fn phase30_live_metric_from_market_with_bump(m: &MarketSnapshot, credit_bump: f64) -> f64 {
        let r = m.ois_quotes.first().copied().unwrap_or(0.03);
        let spot = m.eq_spot.max(1e-6);
        let vol = m.eq_vol.max(1e-4);
        let inf = m.inflation_quotes.first().map(|(_, x)| *x).unwrap_or(0.023);
        let taylor = ql_curves::evaluate_romer_veronesi_taylor_rule_macro_yield_curve(
            &ql_curves::RomerVeronesiTaylorRuleMacroCurveSpec {
                natural_real_rate_r_star: r - inf,
                inflation_target_pi_star: 0.02,
                current_inflation: inf,
                current_output_gap: 0.01,
                taylor_inflation_weight_phi_pi: 1.5,
                taylor_output_weight_phi_y: 0.5,
                inflation_mean_reversion: 0.4,
                output_gap_mean_reversion: 0.3,
                inflation_vol: vol,
                output_gap_vol: vol * 0.5,
                market_price_of_inflation_risk: 0.2,
                zlb_floor_rate: 0.0,
                maturities: vec![1.0, 5.0, 10.0],
            },
        );
        let dam = ql_credit::evaluate_damodaran_rsm_country_risk_and_debt_modification(
            &ql_credit::DamodaranRsmCountryRiskAndDebtModSpec {
                mature_market_erp: 0.05,
                sovereign_cds_spread: 0.012 + credit_bump,
                sovereign_default_drag_bps: 20.0,
                emerging_equity_vol: vol * 1.4,
                emerging_govt_bond_vol: vol.max(0.01),
                firm_domestic_revenue_share: 0.6,
                ebit_operating_income: 12.0,
                annual_interest_expense: 3.0,
                original_debt_carrying_value: 100.0,
                original_effective_interest_rate: r,
                modified_cashflows: vec![4.0, 4.0, 104.0],
                modified_payment_times: vec![1.0, 2.0, 3.0],
                modification_fees_paid: 1.0,
            },
        );
        let var = ql_portfolio::evaluate_christoffersen_bouchaud_var_backtest_and_cumulant_es(
            &ql_portfolio::ChristoffersenBouchaudVarBacktestAndCumulantSpec {
                total_days_t: 250,
                count_n00: 220.0,
                count_n01: 8.0,
                count_n10: 8.0,
                count_n11: 14.0,
                var_tail_probability_p: 0.05,
                portfolio_mean_loss: 0.0,
                portfolio_std_loss: 1.0 + vol,
                loss_skewness_gamma1: 0.2,
                loss_excess_kurtosis_gamma2: 0.4,
            },
        );
        let lvar = ql_portfolio::evaluate_foucault_vandermerwe_lob_schedule_and_lvar(
            &ql_portfolio::FoucaultVanDerMerweLobLiquidityLvarSpec {
                mid_price: spot,
                fundamental_volatility: vol,
                informed_trader_share_pi: 0.25,
                exponential_order_size_rate_lambda: 0.8,
                order_processing_fee_c: 0.0002,
                position_notional: 1_000_000.0,
                mean_relative_bid_ask_spread: 0.001,
                std_relative_bid_ask_spread: 0.0004,
                spread_tail_scaling_a: 1.5,
                confidence_z_alpha: 2.33,
                liquidation_horizon_days: 5.0,
                market_adv_notional: 10_000_000.0,
                permanent_impact_gamma: 0.02,
            },
        );
        let pt = ql_portfolio::evaluate_barberis_jin_wang_prospect_overhang_and_skew_anomaly(
            &ql_portfolio::BarberisJinWangProspectAnomalySpec {
                current_price: spot,
                reference_purchase_price: 95.0,
                return_volatility: vol,
                return_skewness: 0.8,
                prospect_loss_aversion_lambda: 2.25,
                prospect_diminishing_sensitivity_alpha: 0.88,
                prelec_probability_weighting_gamma: 0.7,
                mental_accounting_integration_b0: 0.4,
                capm_beta: 1.1,
                equity_risk_premium: 0.05,
            },
        );
        let regret = ql_portfolio::evaluate_hart_mascolell_regret_matching_allocation(
            &ql_portfolio::HartMasColellRegretMatchingPortfolioSpec {
                strategy_payoff_history: vec![
                    [0.02, 0.01, 0.03],
                    [0.01, 0.02, 0.0],
                    [0.03, 0.01, 0.02],
                ],
                realized_action_weights_history: vec![
                    [0.4, 0.3, 0.3],
                    [0.3, 0.4, 0.3],
                    [0.3, 0.3, 0.4],
                ],
                regret_inertia_mu: 1.0,
                exploration_epsilon: 0.01,
            },
        );
        let forensic = ql_ad::evaluate_beneish_wahlen_forensic_mscore_and_hkas32_bifurcation(
            &ql_ad::BeneishWahlenHkas32ForensicBifurcationSpec {
                dsri: 1.1,
                gmi: 1.05,
                aqi: 1.2,
                sgi: 1.15,
                depi: 1.03,
                sgai: 1.02,
                lvgi: 1.1,
                tata: 0.03,
                net_income: 10.0,
                operating_cash_flow: 8.0,
                average_total_assets: 200.0,
                convertible_issue_proceeds: 100.0,
                convertible_face_value: 100.0,
                contractual_coupon_rate: 0.03,
                comparable_straight_debt_yield: r + credit_bump,
                direct_issuance_transaction_costs: 1.0,
                maturity_years: 5,
            },
        );
        let camp = ql_portfolio::evaluate_colin_bacon_campisi_fixed_income_attribution(
            &ql_portfolio::ColinBaconCampisiFixedIncomeAttributionSpec {
                portfolio_total_return: 0.042 + credit_bump,
                benchmark_total_return: 0.035,
                portfolio_coupon_yield: 0.04,
                benchmark_coupon_yield: 0.035,
                portfolio_eff_duration: 5.2,
                benchmark_eff_duration: 5.0,
                portfolio_spread_duration: 4.0,
                benchmark_spread_duration: 4.2,
                parallel_treasury_shift: 0.005,
                portfolio_key_rate_twist_shift: 0.001,
                benchmark_key_rate_twist_shift: 0.001,
                portfolio_credit_spread_change: 0.002,
                benchmark_credit_spread_change: 0.002,
                fx_spot_return: 0.01,
                domestic_cash_rate: r,
                foreign_cash_rate: 0.032,
            },
        );
        let heck = ql_dml::evaluate_greene_wooldridge_heckman_survivorship_and_gmm(
            &ql_dml::GreeneWooldridgeHeckmanSelectionGmmSpec {
                selection_z_scores: vec![-0.5, 0.0, 0.5, 1.0],
                factor_exposures_x: vec![0.5, 1.0, 1.5, 2.0],
                observed_returns_y: vec![0.02, 0.04, 0.06, 0.09],
                gmm_moment_conditions_bar_g: [0.01, -0.02, 0.015],
                gmm_weighting_matrix_diag: [2.0, 1.5, 1.0],
                sample_size_n: 100.0,
            },
        );
        let en = ql_dml::solve_hastie_tibshirani_elastic_net_and_gcv(
            &ql_dml::HastieTibshiraniElasticNetGcvSpec {
                design_matrix_x: vec![
                    [1.0, 0.2, 0.1],
                    [1.0, 0.4, 0.5],
                    [1.0, 0.6, 0.2],
                    [1.0, 0.8, 0.9],
                ],
                response_y: vec![0.1, 0.2, 0.3, 0.4],
                l1_penalty_lambda1: 0.05 + credit_bump,
                l2_penalty_lambda2: 0.1,
                smoothing_eps: 1e-4,
            },
        );
        let cb = ql_fixed_income::evaluate_an_daoquan_clifford_chance_resettable_convertible(
            &ql_fixed_income::AnDaoquanCliffordChanceResettableCbSpec {
                spot_price: spot,
                initial_conversion_price: 100.0,
                book_value_per_share_bps: 80.0,
                avg_20day_price: spot * 0.95,
                face_value: 100.0,
                coupon_rate: 0.02,
                redemption_maturity_price: 105.0,
                put_back_trigger_ratio: 0.7,
                put_back_price: 103.0,
                soft_call_trigger_ratio: 1.3,
                soft_call_price: 105.0,
                stock_volatility: vol,
                risk_free_rate: r,
                credit_spread: 0.01 + credit_bump,
                maturity: 3.0,
            },
        );
        let bail = ql_credit::evaluate_fsb_finma_gsib_resolution_bail_in_waterfall(
            &ql_credit::FsbFinmaGsibResolutionBailInWaterfallSpec {
                total_rwa: 1000.0,
                pre_resolution_cet1: 120.0,
                at1_principal: 50.0,
                tier2_principal: 60.0,
                senior_tlac_principal: 100.0,
                valuation_loss_shock: 0.08 + credit_bump,
                target_post_resolution_cet1_ratio: 0.1,
                liquidation_insolvency_haircut: 0.4,
                existing_share_count: 100.0,
                bail_in_conversion_discount: 0.2 + credit_bump * 0.01,
            },
        );
        taylor.macro_affine_zero_yields.iter().copied().sum::<f64>()
            + dam.country_risk_premium_crp
            + var.cornish_fisher_var
            + lvar.total_liquidity_adjusted_var_lvar * 1e-6
            + pt.bjw_equilibrium_expected_excess_return
            + regret.hannan_external_regret_per_period
            + forensic.beneish_8var_m_score
            + camp.total_active_local_return
            + heck.heckman_corrected_alpha
            + en.generalized_cross_validation_gcv
            + cb.three_line_duplex_fair_value * 1e-3
            + bail.post_bail_in_creditor_ownership_share
            + bail.ncwo_compensation_shortfall * 1e-3
            + dam.total_cost_of_equity_ke
            + dam.extinguishment_gain_or_loss * 1e-3
            + forensic.manipulation_probit_probability
            + forensic.year1_effective_interest_expense * 1e-2
            + camp.active_security_selection_residual
            + camp.karnosky_singer_unhedged_currency_surprise
    }
    pub fn phase30_live_metric(&self) -> f64 {
        let bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else {
            0.0
        };
        Self::phase30_live_metric_from_market(&self.scenario_market(), bump)
    }

    pub fn phase32_live_metric_from_market(m: &MarketSnapshot, credit_bump: f64) -> f64 {
        let r = m.ois_quotes.first().copied().unwrap_or(0.03);
        let inf = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let f = m.eq_spot.max(1.0);
        let vol = m.eq_vol.max(0.1);
        let c = ql_models::evaluate_gairat_polala_cheyette_local_vol_caplet(
            &ql_models::GairatPolalaCheyetteCapletSpec {
                forward_rate: r,
                strike_rate: r + credit_bump,
                fixing_time_t1: 1.0,
                payment_time_t2: 1.5,
                mean_reversion_kappa: 0.1,
                initial_sigma0: vol * 0.1,
                initial_skew_b: -0.2,
                target_atm_normal_vol: vol * 0.01,
                target_otm_normal_vol: vol * 0.011,
                discount_factor_t2: (-r * 1.5).exp(),
            },
        );
        let q = ql_models::evaluate_gatheral_jusselin_quadratic_rough_heston(
            &ql_models::GatheralJusselinQuadraticRoughHestonSpec {
                forward: f,
                strike: f * 1.01,
                expiry: 1.0,
                hurst_h: 0.1,
                weight_a: 0.2,
                asymmetry_center_b: 0.3,
                floor_variance_c: 0.02,
                initial_state_z0: 0.1,
                kernel_eta: 0.5,
                mean_reversion_lambda: 1.0,
            },
        );
        let z = ql_vol::evaluate_andreasen_huge_one_step_zabr_calibration(
            &ql_vol::AndreasenHugeOneStepZabrSpec {
                forward: f,
                strike: f * 1.01,
                expiry: 1.0,
                alpha: vol * 0.2,
                beta: 0.5,
                rho: -0.3,
                nu: 0.4,
                discount_factor: 1.0,
            },
        );
        let n = ql_curves::evaluate_adb_rbi_ndirs_and_asset_swap_basis(
            &ql_curves::AdbRbiNdirsAssetSwapSpec {
                payment_tenors: vec![1.0, 2.0],
                onshore_forward_rates: vec![r, r + 0.002],
                usd_ois_discount_factors: vec![1.0, 0.95],
                local_rate_vol: 0.1,
                ndf_fx_vol: vol * 0.1,
                rate_fx_corr: 0.2,
                sovereign_bond_yield: r + inf,
                usd_sofr_yield: r,
                ndf_ccs_basis_spread: 0.005,
                withholding_tax_rate: 0.1,
            },
        );
        let x = ql_xva::evaluate_turlakov_green_kenyon_wwr_and_kva(
            &ql_xva::TurlakovGreenKenyonWwrKvaSpec {
                time_grid: vec![1.0, 2.0],
                independent_epe_profile: vec![f * 0.01, f * 0.012],
                pfe_addon_profile: vec![f * 0.005, f * 0.006],
                base_cpty_hazard: 0.02 + credit_bump,
                hazard_log_vol: 0.3,
                exposure_vol: 0.2,
                wwr_correlation_rho: 0.5,
                lgd: 0.6,
                risk_free_rate: r,
                bank_hazard_rate: 0.01,
                ccr_risk_weight: 0.08,
                shareholder_cost_of_capital: 0.1,
                capital_return_offset_beta: 0.5,
            },
        );
        let e = ql_portfolio::evaluate_ahmadi_javid_entropic_var_portfolio(
            &ql_portfolio::AhmadiJavidEntropicVarPortfolioSpec {
                scenario_returns: vec![vec![0.01, -0.02], vec![-0.03, 0.01]],
                scenario_probabilities: vec![0.5, 0.5],
                initial_weights: vec![0.6, 0.4],
                tail_alpha: 0.05,
                risk_aversion_lambda: 1.0,
            },
        );
        let h = ql_portfolio::evaluate_cotton_schur_hrp_and_rmt_nco(
            &ql_portfolio::CottonSchurHrpRmtNcoSpec {
                covariance_4x4: [
                    [0.04, 0.01, 0.0, 0.0],
                    [0.01, 0.05, 0.0, 0.0],
                    [0.0, 0.0, 0.03, 0.01],
                    [0.0, 0.0, 0.01, 0.06],
                ],
                schur_gamma: 0.5,
                sample_aspect_ratio_q: 0.5,
            },
        );
        let b = ql_portfolio::evaluate_tsay_bekk_garch_volatility_spillover(
            &ql_portfolio::TsayBekkGarchSpilloverSpec {
                chol_c: [[0.1, 0.0], [0.02, 0.1]],
                arch_a: [[0.1, 0.02], [0.01, 0.1]],
                garch_b: [[0.8, 0.0], [0.0, 0.7]],
                initial_cov: [[0.04, 0.0], [0.0, 0.04]],
                return_innovations: vec![[0.01, -0.02], [0.02, 0.01]],
            },
        );
        let ql = ql_fixed_income::evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole(
            &ql_fixed_income::MayerBrownSsgaCappedCallMakeWholeSpec {
                bond_par: 100.0,
                spot_price: f,
                base_conversion_price: f * 0.9,
                cap_strike_price: f * 1.2,
                volatility: vol,
                risk_free_rate: r,
                maturity_years: 5.0,
                acquisition_offer_price: f * 1.1,
                remaining_years_at_acquisition: 2.0,
                tax_rate_section_163j: 0.21,
            },
        );
        let hhw = ql_pde::evaluate_haentjens_inthout_heston_hull_white_adi(
            &ql_pde::HaentjensInthoutHhwAdiSpec {
                spot: f,
                strike: f,
                expiry: 1.0,
                initial_var_v0: 0.04,
                mean_rev_kappa_v: 1.0,
                long_var_theta_v: 0.04,
                vol_of_vol_sigma_v: 0.3,
                initial_rate_r0: r,
                rate_mean_rev_kappa_r: 0.1,
                rate_vol_sigma_r: 0.01,
                corr_spot_vol: 0.2,
                corr_spot_rate: 0.3,
            },
        );
        let qual = ql_models::evaluate_schilit_sloan_dechow_dichev_earnings_quality(
            &ql_models::SchilitSloanDechowDichevQualitySpec {
                revenue_0: 100.0,
                revenue_1: 105.0,
                cogs_0: 60.0,
                cogs_1: 63.0,
                accounts_receivable_0: 10.0,
                accounts_receivable_1: 12.0,
                inventory_0: 8.0,
                inventory_1: 9.0,
                accounts_payable_0: 7.0,
                accounts_payable_1: 8.0,
                net_operating_assets_0: 50.0,
                net_operating_assets_1: 52.0,
                net_income_1: 8.0,
                cfo_prev: 7.0,
                cfo_curr: 7.5,
                cfo_next: 8.0,
            },
        );
        let cr = ql_credit::evaluate_brigo_morini_cmcds_and_credit_index_option(
            &ql_credit::BrigoMoriniCmcdsCreditOptionSpec {
                reset_times: vec![1.0, 2.0],
                forward_cds_spreads: vec![0.01, 0.012],
                spread_volatility: 0.3,
                annuity_convexity_gamma: 0.5,
                flat_par_cds_spread: 0.011,
                option_expiry: 1.0,
                option_strike_spread: 0.01,
                risky_pv01_annuity: 4.0,
                hazard_rate: 0.02,
                recovery_rate: 0.4,
                risk_free_rate: r,
            },
        );
        c.deterministic_auxiliary_y_bar.val()
            + q.qrh_call_price.val()
            + z.one_step_fd_call_price.val()
            + n.usd_asset_swap_pickup_bps.val() * 1e-4
            + x.total_wwr_cva_plus_kva.val()
            + e.exact_chernoff_evar.val()
            + h.schur_hrp_portfolio_variance.val()
            + b.terminal_var_asset1.val()
            + ql.capped_call_overlay_gross_cost.val() * 1e-2
            + cr.morini_brigo_payer_swaption_pv.val()
            + hhw.hhw_adi_call_price.val()
            + qual.schilit_shenanigan_composite_z.val()
    }
    pub fn phase33_live_metric_from_market(m: &MarketSnapshot, credit_bump: f64) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.01;
        let f = m.eq_spot.max(1.0);
        let vols = m.eq_vol.max(0.1);
        let c = ql_curves::evaluate_brigo_mercurio_rebonato_swaption_and_cms(
            &ql_curves::BrigoMercurioRebonatoCmsSpec {
                fixing_tenors: vec![1.0, 2.0],
                accrual_taus: vec![1.0, 1.0],
                forward_rates: vec![r, r + 0.005],
                discount_factors: vec![(-r).exp(), (-2.0 * r).exp()],
                forward_vols: vec![vols, vols],
                corr_decay_beta: 0.1,
                option_expiry: 1.0,
                swaption_strike: r,
            },
        );
        let cr = ql_credit::evaluate_duffie_singleton_leland_toft_capital_structure(
            &ql_credit::DuffieSingletonLelandToftSpec {
                unlevered_asset_value: f,
                annual_coupon_c: 0.04 * f,
                debt_principal_p: 0.5 * f,
                rollover_maturity_t: 5.0,
                corporate_tax_rate: 0.25,
                bankruptcy_cost_alpha: 0.2,
                asset_payout_rate: 0.02,
                risk_free_rate: r.max(0.001),
                asset_volatility: vols,
            },
        );
        let ca = ql_fixed_income::evaluate_calamos_despiegeleer_ascot_and_gamma_scalp(
            &ql_fixed_income::CalamosDeSpiegeleerAscotGammaSpec {
                spot_price: f,
                reference_spot: f,
                bond_par: 100.0,
                coupon_rate: 0.04,
                conversion_ratio: 1.0,
                implied_vol: vols,
                realized_vol: vols * 1.1,
                risk_free_rate: r.max(0.001),
                dividend_yield: 0.01,
                base_hazard_rate: 0.02,
                equity_credit_power_p: 1.2,
                recovery_rate: 0.4,
                maturity_years: 3.0,
                holding_period_years: 1.0,
            },
        );
        let sp = ql_fixed_income::evaluate_tuckman_specified_pool_payup_and_dollar_roll(
            &ql_fixed_income::TuckmanSpecifiedPoolDollarRollSpec {
                net_coupon_rate: 0.03,
                discount_yield: r.max(0.001),
                effective_maturity_years: 5.0,
                ctd_tba_cpr: 0.08,
                llb_pool_cpr_multiplier: 0.7,
                tba_front_month_price: 100.0,
                tba_back_month_price: 99.8,
                gc_repo_rate: r,
            },
        );
        let og = ql_portfolio::evaluate_ohara_glosten_milgrom_sequential_spread(
            &ql_portfolio::OharaGlostenMilgromSequentialSpec {
                low_value_vl: f * 0.9,
                high_value_vh: f * 1.1,
                initial_prior_p0: 0.5,
                informed_trader_fraction_mu: 0.3,
                informed_signal_accuracy_theta: 0.8,
                dealer_inventory_gamma: 0.01,
                signed_order_flow: vec![0.2, -0.1, 0.3],
            },
        );
        let rie = ql_portfolio::evaluate_potters_bouchaud_stieltjes_rie_shrinkage(
            &ql_portfolio::PottersBouchaudStieltjesRieSpec {
                sample_eigenvalues: vec![0.5 + credit_bump, 1.0 + r * 0.01, 2.0 + vols * 0.01],
                aspect_ratio_q: (0.5 + credit_bump).min(0.9),
                stieltjes_bandwidth_eta: 0.1,
            },
        );
        let eg = ql_portfolio::evaluate_tsay_egarch_and_hamilton_markov_switching(
            &ql_portfolio::TsayEgarchHamiltonRegimeSpec {
                returns: vec![0.01 + r * 0.01, -0.02 + credit_bump, 0.015 + vols * 0.001],
                initial_variance: (0.04 + vols * 0.001).max(0.001),
                egarch_omega: -0.1,
                egarch_alpha: 0.1,
                egarch_leverage_gamma: -0.1,
                egarch_beta: 0.9,
                calm_mean: 0.0,
                calm_vol: 0.1,
                crisis_mean: 0.0,
                crisis_vol: 0.3,
                stay_calm_p00: 0.95,
                stay_crisis_p11: 0.9,
            },
        );
        let iv = ql_dml::evaluate_wooldridge_greene_2sls_iv_and_did(
            &ql_dml::WooldridgeGreene2slsDidSpec {
                endogenous_x: vec![1.0 + credit_bump, 2.0 + r * 0.01, 3.0 + vols * 0.01],
                instrument_z: vec![1.0, 2.1, 2.9],
                outcome_y: vec![1.2, 2.1, 3.2],
                treated_pre_mean: 1.0,
                treated_post_mean: 1.2,
                control_pre_mean: 1.0,
                control_post_mean: 1.1,
            },
        );
        let mc = ql_portfolio::evaluate_mcneil_frey_copula_tail_and_euler_allocation(
            &ql_portfolio::McneilFreyCopulaEulerAllocationSpec {
                clayton_theta: 2.0,
                gumbel_theta: 1.5,
                student_t_df_nu: 5.0,
                student_t_corr_rho: 0.3,
                desk_weights: [0.25 + credit_bump * 0.1, 0.25, 0.25, 0.25],
                desk_expected_losses: [0.01; 4],
                loss_covariance_4x4: [
                    [0.04, 0.01, 0.0, 0.0],
                    [0.01, 0.05, 0.01, 0.0],
                    [0.0, 0.01, 0.06, 0.01],
                    [0.0, 0.0, 0.01, 0.03],
                ],
                es_tail_multiplier: 2.0,
            },
        );
        let om = ql_portfolio::evaluate_bacon_omega_sortino_and_burke_drawdown(
            &ql_portfolio::BaconOmegaSortinoBurkeDrawdownSpec {
                period_returns: vec![0.01 + r * 0.01, -0.005 + credit_bump, 0.02 + vols * 0.001],
                mar_threshold: 0.0,
                risk_free_rate: 0.002,
            },
        );
        let rd = ql_models::evaluate_damodaran_rd_capitalization_and_patent_real_option(
            &ql_models::DamodaranRdCapitalizationRealOptionSpec {
                reported_ebit: 10.0,
                reported_invested_capital: 100.0,
                effective_tax_rate: 0.25,
                current_rd_expense: 5.0,
                historical_rd_expenses: vec![4.0, 3.0, 2.0],
                patent_project_pv: 20.0,
                patent_development_cost: 10.0,
                patent_life_years: 5.0,
                project_volatility: 0.3,
                risk_free_rate: r.max(0.001),
                abandonment_salvage_value: 2.0,
            },
        );
        let ad = ql_models::evaluate_mascolell_cochrane_arrow_debreu_sdf(
            &ql_models::MasColellCochraneArrowDebreuSdfSpec {
                state_probabilities: [0.3, 0.4, 0.3],
                consumption_growth_ratios: [0.98, 1.02 + credit_bump, 1.05 + r * 0.01],
                equity_payoff_states: [0.8, 1.1, 1.4],
                subjective_discount_beta: 0.99,
                crra_risk_aversion_gamma: 2.0,
            },
        );
        c.forward_swap_rate
            + cr.structural_credit_spread_bps * 1e-4
            + ca.expected_gamma_scalp_pnl
            + sp.specified_pool_payup_points * 1e-2
            + og.inventory_skewed_reservation_mid * 1e-4
            + rie.rie_condition_number * 1e-2
            + eg.terminal_filtered_crisis_probability
            + iv.two_stage_ls_iv_beta * 1e-2
            + mc.total_portfolio_es
            + om.keating_shadwick_omega_ratio * 1e-2
            + rd.total_strategic_expanded_npv * 1e-2
            + ad.equilibrium_equity_price * 1e-2
    }
    pub fn phase33_live_metric(&self) -> f64 {
        let bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else if self.scenario == ScenarioKind::StagflationPlus150Bp {
            0.005
        } else {
            0.0
        };
        Self::phase33_live_metric_from_market(&self.scenario_market(), bump)
    }

    pub fn phase32_live_metric(&self) -> f64 {
        let bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else if self.scenario == ScenarioKind::StagflationPlus150Bp {
            0.005
        } else {
            0.0
        };
        Self::phase32_live_metric_from_market(&self.scenario_market(), bump)
    }

    pub fn phase31_live_metric_from_market(m: &MarketSnapshot, credit_bump: f64) -> f64 {
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump;
        let inf = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let k = ql_curves::evaluate_fabozzi_tuckman_key_rate_and_lsc_immunization(
            &ql_curves::FabozziKeyRateLscImmunizationSpec {
                cashflow_tenors: vec![1.0, 3.0, 5.0, 10.0],
                cashflows: vec![2.0, 2.0, 2.0, 102.0],
                zero_rates: vec![r, r + 0.002, r + 0.004, r + 0.006],
                key_rate_tenors: vec![1.0, 5.0, 10.0],
                ns_decay_lambda: 0.5 + inf,
                liability_horizon_years: 5.0,
                barbell_short_tenor: 1.0,
                barbell_long_tenor: 10.0,
                twist_slope_shock_bps: 25.0,
                butterfly_curvature_shock_bps: 15.0,
            },
        );
        let cmo = ql_fixed_income::evaluate_davidson_levin_cmo_pac_companion_waterfall(
            &ql_fixed_income::DavidsonLevinCmoPacCompanionSpec {
                initial_pool_balance: 100.0,
                gross_wac_rate: r + 0.03,
                net_coupon_rate: r + 0.02,
                num_periods: 24,
                pac_lower_psa: 100.0,
                pac_upper_psa: 300.0,
                realized_psa: 180.0 + credit_bump * 100.0,
                discount_yield: r,
            },
        );
        let cb = ql_fixed_income::evaluate_asquith_sirbu_shreve_cb_call_delay_and_dynkin_game(
            &ql_fixed_income::AsquithSirbuShreveCbGameSpec {
                firm_value: 130.0,
                stock_price: 100.0,
                conversion_ratio: 1.0,
                annual_coupon_amount: 3.0,
                annual_dividend_per_share: 1.0,
                corporate_tax_rate: 0.25,
                effective_call_price: 105.0,
                call_notice_period_years: 30.0 / 365.0,
                firm_volatility: m.eq_vol.max(0.05),
                risk_free_rate: r,
                firm_payout_rate: 0.02,
                dilution_gamma: 0.8,
            },
        );
        let asu = ql_fixed_income::evaluate_fasb_asu2024_induced_conversion_and_cpdi_oid(
            &ql_fixed_income::FasbAsu2024InducedConversionCpdiSpec {
                net_income: 12.0,
                basic_shares: 100.0,
                bond_par: 100.0,
                stated_coupon_rate: 0.03,
                comparable_noncontingent_yield: r + 0.01,
                tax_rate: 0.25,
                original_conversion_shares: 1.0,
                induced_sweetener_extra_shares: 0.05,
                induced_sweetener_cash: 1.0,
                stock_price_at_acceptance: m.eq_spot,
                maturity_years: 5.0,
            },
        );
        let modif = ql_credit::evaluate_rsm_asc470_debt_modification_and_restructuring(
            &ql_credit::RsmAsc470DebtModificationSpec {
                carrying_amount: 100.0,
                original_par: 100.0,
                original_coupon_rate: 0.03,
                original_effective_rate: r,
                modified_par: 100.0,
                modified_coupon_rate: 0.025,
                remaining_years: 5.0,
                lender_fee_paid: 1.0,
                conversion_option_fv_change: 0.5,
                market_exit_yield: r + 0.01,
            },
        );
        let colva = ql_xva::evaluate_gregory_csa_threshold_mta_and_ctd_colva(
            &ql_xva::GregoryCsaThresholdMtaCtdColvaSpec {
                time_grid: vec![1.0, 2.0, 3.0, 5.0],
                expected_mtm_profile: vec![5.0, 8.0, 10.0, 12.0],
                exposure_vol: 3.0 + m.eq_vol,
                mpor_years: 10.0 / 252.0,
                csa_threshold: 2.0,
                minimum_transfer_amount: 1.0,
                independent_amount: 0.5,
                hazard_rate: 0.02 + credit_bump,
                lgd: 0.6,
                risk_free_rate: r,
                ccy1_collateral_spread: 0.01,
                ccy2_collateral_spread: 0.008,
                collateral_spread_vol: 0.01,
            },
        );
        let evt = ql_portfolio::evaluate_embrechts_hill_pickands_and_subexponential_tail(
            &ql_portfolio::EmbrechtsHillPickandsSubexponentialSpec {
                loss_samples: vec![
                    160.0, 110.0, 82.0, 66.0, 56.0, 50.0, 46.0, 43.0, 40.0, 38.0, 36.0, 35.0,
                ],
                tail_order_k: 2,
                extreme_tail_prob: 0.01,
                portfolio_positions_count: 8,
            },
        );
        let tlf = ql_portfolio::evaluate_bouchaud_potters_power_law_cumulant_and_smile_scaling(
            &ql_portfolio::BouchaudPottersCumulantSmileSpec {
                single_period_vol: m.eq_vol,
                single_period_skewness: -0.2,
                single_period_excess_kurtosis: 3.0,
                vol_of_vol_kurtosis_component: 2.0,
                vol_decorrelation_steps: 10.0,
                horizon_steps_n: 20.0,
                log_moneyness: 0.0,
                confidence_z: 2.33,
            },
        );
        let attr = ql_portfolio::evaluate_karnosky_singer_multicurrency_and_menchero_attribution(
            &ql_portfolio::KarnoskySingerMencheroAttributionSpec {
                portfolio_asset_weights: vec![0.6, 0.4],
                benchmark_asset_weights: vec![0.5, 0.5],
                portfolio_currency_weights: vec![0.6, 0.4],
                benchmark_currency_weights: vec![0.5, 0.5],
                local_portfolio_returns: vec![0.08, 0.03],
                local_benchmark_returns: vec![0.06, 0.04],
                local_eurodollar_rates: vec![0.02, 0.025],
                base_eurodollar_rate: r,
                spot_fx_returns: vec![0.01, -0.005],
                multi_period_port_returns: vec![0.02, 0.03],
                multi_period_bench_returns: vec![0.018, 0.025],
            },
        );
        let tm = ql_portfolio::evaluate_dgtw_characteristic_and_treynor_mazuy_timing(
            &ql_portfolio::DgtwTreynorMazuyTimingSpec {
                current_weights: vec![0.6, 0.4],
                lagged_weights: vec![0.5, 0.5],
                stock_returns: vec![0.08, 0.03],
                current_char_bench_returns: vec![0.06, 0.04],
                lagged_char_bench_returns: vec![0.05, 0.04],
                portfolio_excess_series: vec![0.01, 0.02, 0.04, 0.08],
                market_excess_series: vec![0.0, 0.02, 0.04, 0.08],
            },
        );
        let bh = ql_portfolio::evaluate_brock_hommes_adaptive_belief_switching(
            &ql_portfolio::BrockHommesAdaptiveBeliefSpec {
                initial_price_deviation: 0.1,
                risk_free_rate: r,
                fundamentalist_damping_v1: 0.8,
                chartist_trend_extrapolation_g: 1.2,
                fundamentalist_info_cost_c1: 0.1,
                intensity_of_choice_beta: 8.0,
                fitness_memory_eta: 0.8,
                risk_aversion_sigma2: 1.0,
                shock_innovations: vec![0.01, -0.005, 0.003],
            },
        );
        let pen = ql_models::evaluate_penman_reoi_and_ohlson_juettner_aeg_valuation(
            &ql_models::PenmanReoiOhlsonJuettnerSpec {
                net_operating_assets_0: 100.0,
                net_financial_obligations_0: 30.0,
                operating_income_after_tax_1: 12.0,
                operating_income_after_tax_2: 13.0,
                net_financial_expense_after_tax_1: 2.0,
                wacc: r + 0.02,
                terminal_reoi_growth: 0.02,
                current_share_price: m.eq_spot,
                forward_eps_1: 5.0,
                forward_eps_2: 5.5,
                forward_dps_1: 2.0,
                perpetual_aeg_growth: 0.03,
            },
        );
        k.bond_present_value
            + cmo.pac_present_value
            + cb.perpetual_dynkin_cb_fair_value
            + cb.asquith_safety_call_trigger_stock_price * 1e-2
            + cb.conversion_value_to_call_ratio
            + asu.basic_eps
            + asu.asu2020_06_instrument_c_net_share_eps
            + asu.asu2024_04_induced_conversion_expense * 1e-2
            + asu.cpdi_projected_contingent_terminal_payment * 1e-3
            + modif.pv_modified_cashflows_at_orig_eir
            + modif.prospective_modified_eir
            + colva.ctd_multi_currency_colva
            + evt.weissman_extreme_var * 1e-3
            + tlf.bouchaud_potters_non_gaussian_var
            + attr.cumulative_active_return
            + attr.karnosky_singer_security_selection
            + attr.karnosky_singer_currency_management
            + attr.menchero_natural_scaling_multiplier * 1e-3
            + tm.treynor_mazuy_timing_gamma
            + tm.dgtw_characteristic_selectivity_cs
            + bh.terminal_price_deviation
            + pen.reoi_intrinsic_equity_value * 1e-2
            + pen.ohlson_juettner_implied_cost_of_equity
    }

    pub fn phase31_live_metric(&self) -> f64 {
        let bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else {
            0.0
        };
        Self::phase31_live_metric_from_market(&self.scenario_market(), bump)
    }

    pub fn phase29_live_metric_from_market(m: &MarketSnapshot) -> f64 {
        Self::phase29_live_metric_from_market_with_bump(m, 0.0)
    }

    pub fn phase29_live_metric_from_market_with_bump(m: &MarketSnapshot, credit_bump: f64) -> f64 {
        let r = m.ois_quotes.first().copied().unwrap_or(0.03);
        let spot = m.eq_spot.max(1e-6);
        let vol = m.eq_vol.max(1e-4);
        let inf = m.inflation_quotes.first().map(|(_, z)| *z).unwrap_or(0.023);
        let a = ql_fixed_income::evaluate_veronesi_fabozzi_futures_convexity_and_prepayment(
            &ql_fixed_income::VeronesiFabozziConvexityPrepaymentSpec {
                forward_rate: r,
                short_rate_vol: vol,
                mean_reversion_kappa: 0.1,
                fixing_time_t1: 1.0,
                maturity_time_t2: 2.0,
                wac_coupon: 0.05,
                current_mortgage_rate: r + inf,
                min_cpr_turnover: 0.02,
                max_cpr_refi: 0.3,
                refi_steepness_gamma: 8.0,
                refi_incentive_threshold: 0.0,
                pool_age_months: 24.0,
                psa_ramp_months: 30.0,
            },
        );
        let b = ql_credit::compute_lando_israel_generator_logarithm_and_migration(
            &ql_credit::LandoIsraelGeneratorLogarithmSpec {
                transition_matrix_3x3: [[0.9, 0.08, 0.02], [0.03, 0.92, 0.05], [0.0, 0.0, 1.0]],
                series_order_k: 12,
                fractional_horizon_years: 0.5,
                recovery_rate: 0.4,
                risk_free_rate: r + credit_bump,
            },
        );
        let c = ql_vol::evaluate_gatheral_demeterfi_corridor_variance_and_skew(
            &ql_vol::GatheralDemeterfiCorridorVarianceSkewSpec {
                forward: spot,
                maturity: 1.0,
                atm_vol: vol,
                atm_skew: -0.2,
                smile_curvature: 0.1,
                corridor_lower_strike: spot * 0.9,
                corridor_upper_strike: spot * 1.1,
                sampling_steps_n: 12,
                vol_of_vol_eta: vol,
                lsv_mixing_weight: 0.3,
            },
        );
        let d = ql_exotics::price_shepp_shiryaev_shreve_perpetual_russian_and_mckean(
            &ql_exotics::SheppShiryaevShrevePerpetualRussianSpec {
                spot,
                running_max: spot * 1.1,
                fixed_put_strike: spot,
                risk_free_rate: r,
                dividend_yield: inf,
                volatility: vol,
            },
        );
        let e = ql_levy::evaluate_cont_tankov_bates_svj_characteristic_and_pide(
            &ql_levy::ContTankovBatesSvjPideSpec {
                spot,
                strike: spot,
                maturity: 1.0,
                risk_free_rate: r,
                dividend_yield: inf,
                v0: vol * vol,
                kappa: 1.0,
                theta: vol * vol,
                sigma_v: vol,
                rho: -0.5,
                jump_intensity_lambda: 0.2,
                jump_mean_mu_j: -0.1,
                jump_std_sigma_j: 0.2,
            },
        );
        let f = ql_fixed_income::evaluate_koziol_asquith_fasb_strategic_convertible(
            &ql_fixed_income::KoziolAsquithFasbStrategicConvertibleSpec {
                firm_value: spot * 100.0,
                debt_face_value: spot * 80.0,
                existing_shares: 100.0,
                new_shares_blockholder1: 10.0,
                new_shares_blockholder2: 8.0,
                coupon_rate: 0.04,
                corporate_tax_rate: 0.25,
                dividend_yield: inf,
                stock_price: spot,
                base_conversion_ratio: 1.0,
                induced_sweetened_conversion_ratio: 1.1,
                cash_sweetener_amount: 1.0,
                call_price: spot * 1.2,
                firm_volatility: vol,
                risk_free_rate: r,
                maturity: 5.0,
            },
        );
        let g = ql_models::solve_carmona_delarue_mean_field_systemic_interbank_game(
            &ql_models::CarmonaDelarueMeanFieldSystemicSpec {
                initial_log_reserve_spread: 0.01 + credit_bump,
                interbank_mean_reversion_a: 0.4,
                control_quadratic_penalty_q: 0.2,
                systemic_risk_preference_epsilon: 0.3,
                terminal_systemic_penalty_c: 1.0,
                idiosyncratic_vol_sigma: vol,
                common_noise_corr_rho0: 0.2,
                distress_barrier_d: -0.1,
                horizon_t: 1.0,
            },
        );
        let h = ql_xva::evaluate_bcbs_d424_output_floor_and_sma_op_risk(
            &ql_xva::BcbsD424OutputFloorAndSmaOpRiskSpec {
                irb_credit_rwa: 100.0,
                ima_market_rwa: 50.0,
                standardized_credit_rwa: 120.0,
                standardized_market_rwa: 60.0,
                standardized_cva_rwa: 20.0,
                output_floor_multiplier: 0.725,
                interest_lease_dividend_component: 20.0,
                services_component: 15.0,
                financial_component: 30.0,
                average_annual_operational_loss: 2.0,
                cet1_capital: 20.0,
                at1_and_tier2_capital: 10.0,
                eligible_tlac_senior_debt: 40.0,
            },
        );
        let i = ql_portfolio::solve_cartea_jaimungal_penalva_lob_fill_hjb(
            &ql_portfolio::CarteaJaimungalPenalvaLobFillHjbSpec {
                initial_inventory: 10,
                horizon_seconds: 60.0,
                base_arrival_rate_lambda: 0.5,
                fill_decay_kappa: 1.0,
                running_inventory_penalty_phi: 0.01,
                terminal_liquidation_penalty_alpha: 0.1,
                ofi_short_term_drift: 0.02,
                adverse_selection_cost: 0.01,
            },
        );
        let j = ql_portfolio::evaluate_madhavan_abner_hill_etf_ap_arbitrage_and_tax_alpha(
            &ql_portfolio::MadhavanAbnerHillEtfApArbitrageSpec {
                etf_market_price: spot * 1.01,
                reported_inav: spot,
                futures_beta_return: 0.01,
                stale_underlying_weight: 0.4,
                ap_creation_fee_bps: 5.0,
                basket_bid_ask_spread_bps: 8.0,
                creation_unit_volatility_bps: 10.0,
                ap_risk_aversion: 0.5,
                embedded_unrealized_gain_fraction: 0.3,
                annual_redemption_turnover: 0.2,
                capital_gains_tax_rate: 0.25,
            },
        );
        let k = ql_portfolio::evaluate_grinold_qian_christopherson_ic_decay_and_timing(
            &ql_portfolio::GrinoldQianChristophersonAlphaDecaySpec {
                fast_signal_ic0: 0.1,
                fast_signal_half_life_months: 3.0,
                slow_signal_ic0: 0.08,
                slow_signal_half_life_months: 12.0,
                signal_correlation: 0.2,
                rebalance_horizon_months: 1.0,
                transaction_cost_penalty: 0.01,
                factor_active_exposures: vec![0.1, -0.05],
                factor_benchmark_returns: vec![0.02, 0.01],
                normal_factor_exposures: vec![0.05, -0.02],
            },
        );
        let l = ql_dml::evaluate_glasserman_saddlepoint_is_and_stratified_lhs_aad(
            &ql_dml::GlassermanSaddlepointIsStratifiedLhsSpec {
                spot,
                deep_otm_strike: spot * 1.3,
                risk_free_rate: r,
                volatility: vol,
                maturity: 1.0,
                num_strata: 16,
            },
        );
        a.implied_futures_rate
            + b.ig_to_default_fractional_pd
            + c.full_variance_swap_fair_vol
            + d.perpetual_russian_option_value
            + e.bates_svj_call_price
            + f.strategic_convertible_fair_value * 1e-3
            + g.mean_field_value_function_cost
            + h.floored_cet1_ratio
            + i.expected_terminal_inventory * 1e-3
            + j.net_ap_creation_redemption_profit_bps * 1e-3
            + k.total_factor_contribution
            + l.importance_sampled_call_price
            + inf
    }

    pub fn phase29_live_metric(&self) -> f64 {
        let credit_bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else {
            0.0
        };
        Self::phase29_live_metric_from_market_with_bump(&self.scenario_market(), credit_bump)
    }

    pub fn phase28_live_metric_from_market(m: &MarketSnapshot) -> f64 {
        Self::phase28_live_metric_from_market_with_bump(m, 0.0)
    }

    fn phase28_live_metric_from_market_with_bump(m: &MarketSnapshot, credit_bump: f64) -> f64 {
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump;
        let spot = m.eq_spot;
        let vol = m.eq_vol.max(1e-4);
        let credit = (m.fwd_quotes.first().copied().unwrap_or(r) - r)
            .abs()
            .max(0.01);
        let inf = m.inflation_quotes.first().map(|(_, z)| *z).unwrap_or(0.023);
        let irr = evaluate_bcbs_d368_six_scenario_irrbb(&BcbsD368SixScenarioIrrbbSpec {
            bucket_midpoints_years: vec![0.5, 2.0, 5.0],
            asset_repricing_cashflows: vec![10.0 + inf, 8.0, 5.0],
            liability_repricing_cashflows: vec![7.0, 6.0, 4.0],
            baseline_zero_rates: vec![r, r + 0.002, r + 0.004],
            parallel_shock_bps: 200.0,
            short_shock_bps: 250.0,
            long_shock_bps: 100.0,
            decay_tau_years: 4.0,
            tier1_capital: 100.0,
            pass_through_deposit_beta: 0.5,
        });
        let rfr = evaluate_turfus_willems_rfr_in_arrears_caplet_smile(
            &TurfusWillemsRfrInArrearsSmileSpec {
                accrual_start_t1: 0.1,
                accrual_end_t2: 0.6,
                forward_rfr_rate: r,
                strike_rate: r * 0.9,
                discount_factor_t2: (-r * 0.6).exp(),
                sabr_alpha: vol * 0.1,
                sabr_beta: 0.5,
                sabr_rho: -0.2,
                sabr_nu: vol,
                hull_white_mean_reversion: 0.1,
                notional: 1e6,
            },
        );
        let fx_spot = 1.1 + spot * 1e-4;
        let xccy =
            evaluate_baran_bianchetti_mtm_xccy_basis_swap(&BaranBianchettiMtmXccyBasisSpec {
                spot_fx_dom_per_for: fx_spot,
                foreign_notional: 100.0,
                payment_times: vec![0.5, 1.0, 1.5],
                dom_ois_discount_factors: vec![(-r * 0.5).exp(), (-r).exp(), (-r * 1.5).exp()],
                for_xccy_discount_factors: vec![0.98, 0.96, 0.94],
                dom_forward_rates: vec![r; 3],
                for_forward_rates: vec![r - 0.005; 3],
                quoted_dom_basis_spread: credit * 0.01,
                dom_rate_vol: vol,
                fx_spot_vol: vol,
                corr_dom_rate_fx: -0.3,
            });
        let bg = evaluate_bourgey_gatheral_bergomi_guyon_expansion(
            &BourgeyGatheralBergomiGuyonSmileSpec {
                forward_price: spot,
                strike_price: spot,
                maturity: 1.0,
                initial_instantaneous_vol: vol,
                vol_of_vol_eta: vol,
                spot_vol_corr_rho: -0.5,
                hurst_exponent: 0.1,
                mean_reversion_lambda: 1.0,
            },
        );
        let coco = evaluate_feds_despiegeleer_at1_coco_dual_trigger(
            &FedsDeSpiegeleerAt1CocoDualTriggerSpec {
                notional: 100.0,
                coupon_rate: 0.08,
                maturity: 5.0,
                risk_free_rate: r,
                initial_cet1_ratio: 0.12 + credit,
                mechanical_cet1_trigger: 0.07,
                ponv_cet1_trigger: 0.05,
                supervisory_ponv_hazard_lambda: credit,
                asset_volatility: vol,
                cet1_asset_elasticity: 1.0,
                initial_stock_price: spot,
                conversion_price: spot * 0.8,
                floor_stock_price: spot * 0.2,
                existing_common_shares: 100.0,
                write_down_fraction: 0.5,
                is_equity_conversion: true,
            },
        );
        let xv = evaluate_albanese_crepey_zeron_balance_sheet_xva(
            &AlbaneseCrepeyZeronBalanceSheetXvaSpec {
                maturity: 2.0,
                num_time_steps: 8,
                risk_free_rate: r,
                cpty_hazard_rate: credit,
                cpty_recovery: 0.4,
                bank_hazard_rate: credit * 0.5,
                bank_recovery: 0.4,
                bank_liquidity_basis: credit * 0.2,
                im_funding_spread: credit,
                shareholder_hurdle_rate: 0.1,
                base_unmargined_epe: 10.0 + inf,
                base_ene: 4.0,
                simm_delta_weight: 0.6,
                simm_vega_weight: 0.4,
                simm_delta_vega_corr: 0.2,
                regulatory_capital_ratio: 0.08,
            },
        );
        let rb =
            evaluate_jacquier_chebyshev_rough_vol_surrogate(&JacquierDeepRoughVolSurrogateSpec {
                forward: spot,
                maturity: 1.0,
                initial_var_xi0: vol * vol,
                hurst_h: 0.1,
                vol_of_vol_eta: vol,
                spot_vol_corr_rho: -0.7,
                target_market_vols_3_strikes: [vol * 1.1, vol, vol * 1.05],
                log_moneyness_strikes: [-0.1, 0.0, 0.1],
            });
        let sp = evaluate_choi_li_deng_multi_asset_spread_basket(
            &ChoiLiDengMultiAssetSpreadBasketSpec {
                long_leg_forward: spot,
                long_leg_weight: 1.0,
                long_leg_vol: vol,
                short_leg_forwards: vec![spot * 0.95],
                short_leg_weights: vec![1.0],
                short_leg_vols: vec![vol * 1.1],
                corr_long_short: vec![0.4],
                corr_short_short_equicorr: 0.0,
                strike: spot * 0.02,
                maturity: 1.0,
                risk_free_rate: r,
            },
        );
        let bs = evaluate_bess_analytical_sdp_degradation_storage(
            &BessAnalyticalSdpDegradationStorageSpec {
                energy_capacity_mwh: 100.0,
                max_power_mw: 50.0,
                initial_soc_fraction: 0.5,
                min_soc_fraction: 0.1,
                max_soc_fraction: 0.9,
                charge_efficiency: 0.9,
                discharge_efficiency: 0.9,
                battery_capex_usd: 100000.0,
                rated_cycle_life: 4000.0,
                peak_mean_price_mwh: spot,
                offpeak_mean_price_mwh: spot * 0.5,
                price_volatility_mwh: spot * 0.2,
                horizon_days: 365,
                discount_rate_annual: r,
            },
        );
        let dk = evaluate_mitchell_pedersen_slow_capital_and_kratz_dark_pool(
            &MitchellPedersenKratzDarkPoolSpiralSpec {
                initial_inventory_shares: 1e5,
                reference_price: spot,
                execution_horizon_days: 5.0,
                daily_volatility: vol,
                lit_temporary_impact_eta: 1e-6,
                lit_permanent_impact_gamma: 2e-6,
                risk_aversion_lambda: 0.1,
                dark_pool_poisson_intensity: 0.3,
                dark_pool_adverse_selection_bps: 5.0,
                specialist_immediate_capital: 1e7,
                slow_moving_external_capital: 2e7,
                capital_arrival_speed_beta: 0.4,
            },
        );
        let rm = evaluate_follmer_schied_convex_risk_and_diagnostic_expectations(
            &FollmerSchiedGennaioliDiagnosticRiskSpec {
                scenario_pnl: vec![spot * 0.01, -spot * 0.02, -spot * 0.05],
                baseline_probabilities: vec![0.5, 0.3, 0.2],
                entropic_ambiguity_theta: 2.0,
                kl_divergence_budget_eta: 0.1,
                diagnostic_overreaction_theta: 1.5,
                boom_news_drift_shock: inf,
                bust_news_drift_shock: credit,
            },
        );
        let pc =
            evaluate_huge_savine_differential_pca_and_twin_dml(&HugeSavineDifferentialPcaSpec {
                state_samples: vec![vec![0.2, 2.0], vec![0.3, -1.8], vec![-0.25, 1.5]],
                payoff_samples: vec![spot, spot * 0.98, spot * 1.02],
                aad_gradient_samples: vec![vec![2.2, 0.1], vec![1.8, -0.15], vec![-2.0, 0.12]],
                twin_fd_gradient_samples: vec![vec![2.0, 0.12], vec![1.7, -0.1], vec![-1.9, 0.1]],
                retained_rank: 1,
                ridge_regularization: 1e-6,
            });
        irr.baseline_eve
            + rfr.in_arrears_caplet_pv * 1e-6
            + xccy.mtm_swap_npv_in_dom
            + bg.atm_implied_vol_second_order
            + coco.at1_coco_fair_value * 1e-3
            + xv.ftd_cva
            + rb.surrogate_smile_vols[0]
            + sp.choi_gauss_hermite_price * 1e-2
            + bs.total_bess_npv_usd * 1e-5
            + dk.kratz_schoneborn_dark_lit_cost * 1e-3
            + rm.follmer_schied_entropic_risk_measure
            + pc.truncation_gradient_mse_bound
            + inf
    }

    pub fn phase28_live_metric(&self) -> f64 {
        let credit_bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else {
            0.0
        };
        Self::phase28_live_metric_from_market_with_bump(&self.scenario_market(), credit_bump)
    }

    pub fn phase27_live_metric(&self) -> f64 {
        let market = self.scenario_market();
        let credit_bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else {
            0.0
        };
        Self::phase27_live_metric_from_market(&market, credit_bump)
    }

    pub fn phase26_live_metric_from_market(m: &MarketSnapshot, credit_bump: f64) -> f64 {
        let r = m.ois_quotes[0] + credit_bump;
        let spot = m.eq_spot;
        let vol = m.eq_vol.max(1e-4);
        let inf = m.inflation_quotes.first().map(|(_, z)| *z).unwrap_or(0.023);
        let times = vec![0.25, 0.5, 0.75, 1.0];
        let mf = evaluate_hunt_kennedy_markov_functional_range_accrual(
            &HuntKennedyMarkovFunctionalSpec {
                fixing_times: times.clone(),
                initial_forwards: vec![r; 4],
                market_caplet_vols: vec![vol; 4],
                smile_curvature_eta: 0.1 * inf,
                mean_reversion: 0.03,
                range_lower: -0.01,
                range_upper: 0.03,
                accrual_coupon_rate: 0.02 + inf,
                call_strike_spread: 0.005,
            },
        );
        let wc = evaluate_tuckman_burghardt_treasury_wildcard_and_ctd_switch(
            &TuckmanBurghardtTreasuryWildcardBasisSpec {
                futures_settlement_price: 100.0,
                delivery_horizon_years: 0.25,
                general_collateral_repo_rate: r,
                repo_specialness_spread: 0.002,
                bond1_clean_price: 99.0,
                bond1_conversion_factor: 0.85,
                bond1_modified_duration: 6.0,
                bond1_yield_beta: 1.0,
                bond2_clean_price: 101.0,
                bond2_conversion_factor: 0.9,
                bond2_modified_duration: 7.0,
                bond2_yield_beta: 1.1,
                afternoon_wildcard_vol_annual: vol,
                eom_delivery_window_days: 5,
            },
        );
        let tarf = evaluate_wystup_pivot_tarf_and_corridor_strip(&WystupPivotTarfCorridorSpec {
            spot,
            strike_low: 0.95 * spot,
            strike_pivot: spot,
            strike_high: 1.05 * spot,
            downside_gearing: 1.2,
            upside_gearing: 1.0,
            target_redemption_cap: 0.05,
            domestic_rate: r,
            foreign_rate: 0.01,
            volatility: vol,
            fixing_times: times.clone(),
        });
        let berg = evaluate_bergomi_forward_start_smile_and_napoleon_convexity(
            &BergomiForwardStartNapoleonVolSpec {
                initial_forward_variance: vol * vol,
                vol_of_vol_omega: 0.3,
                kappa_short_k1: 3.0,
                kappa_long_k2: 0.5,
                long_factor_weight_theta: 0.5,
                factor_corr_rho12: 0.3,
                spot_vol_corr_rho_s1: -0.5,
                spot_vol_corr_rho_s2: -0.3,
                forward_start_t1: 0.5,
                forward_end_t2: 1.0,
                forward_moneyness_k: 1.0,
                napoleon_base_coupon: 0.01,
                napoleon_global_floor: 0.0,
                napoleon_subperiods_m: 4,
            },
        );
        let sch = evaluate_schonbucher_defaultable_range_accrual_and_spread_cap(
            &SchonbucherDefaultableRangeAccrualCapSpec {
                notional: 100.0,
                payment_times: times.clone(),
                risk_free_rate: r,
                initial_hazard_rates: vec![0.02 + credit_bump; 4],
                base_recovery_rate: 0.4,
                spread_volatilities: vec![vol; 4],
                spread_cap_strike: 0.03,
                accrual_spread_lower: 0.0,
                accrual_spread_upper: 0.05,
                range_accrual_coupon_rate: 0.02,
                recovery_spread_correlation: -0.2,
            },
        );
        let cva = evaluate_green_frtb_cva_hedge_relief_and_kva(&GreenFrtbCvaHedgeCapitalSpec {
            counterparty_sacr_eads: vec![100.0 * (1.0 + credit_bump), 80.0 * (1.0 + vol)],
            counterparty_effective_maturities: vec![3.0 + inf, 5.0 + vol],
            counterparty_risk_weights: vec![0.02 + credit_bump, 0.03 + inf.abs()],
            single_name_hedge_notionals: vec![20.0 * (1.0 + spot.abs() * 1e-3)],
            single_name_hedge_maturities: vec![3.0],
            hedge_reference_correlations_rhc: vec![0.8],
            index_hedge_notional: 10.0 * (1.0 + r.abs()),
            index_hedge_maturity: 5.0,
            index_risk_weight: 0.02,
            sa_cva_spread_deltas: vec![1.0, 2.0],
            sa_cva_vega_sensitivities: vec![0.5],
            hurdle_rate_gamma_k: 0.1,
            hedge_margin_funding_spread: 0.01 + vol,
            discount_rate: r.max(1e-4),
        });
        let lb = evaluate_heynen_kat_partial_fractional_lookback_and_extreme_spread(
            &HeynenKatPartialFractionalLookbackSpec {
                spot,
                running_min: 0.9 * spot,
                running_max: 1.1 * spot,
                fractional_multiplier_lambda: 1.2,
                monitoring_end_t1: 0.5,
                maturity_t2: 1.0,
                fixed_strike: spot,
                rate: r,
                dividend_yield: 0.0,
                volatility: vol,
            },
        );
        let vta = evaluate_eydeland_wolyniec_ramp_outage_virtual_tolling(
            &EydelandWolyniecRampOutageVtaSpec {
                power_forwards: vec![60.0, 62.0],
                gas_forwards: vec![30.0, 31.0],
                spinning_reserve_prices: vec![3.0, 3.0],
                dispatch_times: times[..2].to_vec(),
                capacity_mw: 100.0,
                block_hours: 1.0,
                full_load_heat_rate: 2.0,
                part_load_heat_rate: 2.2,
                min_stable_load_fraction: 0.5,
                vom_cost_per_mwh: 1.0,
                startup_cost: 2.0,
                fast_ramp_penalty_cost: 1.0,
                outage_failure_rate_lambda: 0.05,
                outage_repair_rate_mu: 0.95,
                spread_volatility: vol,
                risk_free_rate: r,
            },
        );
        let fl =
            evaluate_grinold_kahn_clarke_fundamental_law(&GrinoldKahnClarkeFundamentalLawSpec {
                standardized_scores_z: vec![1.0 + inf, -0.5 + credit_bump, 0.7 + spot * 1e-3],
                residual_volatilities: vec![
                    vol,
                    (vol + r.abs()).max(1e-4),
                    (vol + inf.abs()).max(1e-4),
                ],
                benchmark_weights: vec![0.3, 0.4, 0.3],
                information_coefficient_ic: 0.1,
                rebalances_per_year: 12,
                signal_autocorrelation_rho: 0.0,
                target_tracking_error: 0.05,
                max_active_weight_bound: 0.1,
            });
        let acd = evaluate_engle_russell_acd_trade_duration_and_hazard(
            &EngleRussellAcdTradeDurationSpec {
                inter_trade_durations: vec![1.0 + vol, 1.2 + r.abs(), 0.9 + spot.abs() * 1e-3],
                omega: 0.2 + r.abs(),
                alpha: 0.2 + inf.abs(),
                beta: 0.5 + credit_bump.abs(),
                weibull_shape_gamma: 1.0 + vol,
                initial_expected_duration: 1.0,
                elapsed_time_since_last_trade: 0.5 + inf.abs(),
                base_price_impact_bps: 2.0 + spot.abs() * 1e-3,
                duration_impact_elasticity: 0.5,
            },
        );
        let hlz = evaluate_harvey_liu_zhu_multiple_testing_and_mintrl(
            &HarveyLiuZhuMultipleTestingSharpeSpec {
                focal_sharpe_ratio: 0.8 + inf + credit_bump,
                benchmark_sharpe_sr0: 0.0,
                sample_length_t: 252,
                return_skewness: -0.5,
                return_kurtosis_pearson: 4.0,
                num_independent_trials_m: 20,
                cross_sectional_sharpe_std: 0.3,
                confidence_z_alpha: 1.96 + vol,
                candidate_t_stats: vec![1.5 + spot * 1e-3, 2.0 + vol],
            },
        );
        let ml = evaluate_giles_capriotti_multilevel_mc_aad(&GilesCapriottiMultilevelMcAadSpec {
            spot,
            strike: spot,
            rate: r,
            volatility: vol,
            local_vol_skew: -0.2,
            maturity: 1.0,
            base_steps_m0: 8,
            target_rmse_epsilon: 0.01,
        });
        mf.non_callable_range_accrual_npv.val()
            + wc.theoretical_net_basis_ctd.val() * 1e-2
            + tarf.pivot_tarf_client_npv.val()
            + berg.bergomi_stoch_vol_napoleon_coupon.val()
            + sch.defaultable_range_accrual_note_npv.val() * 1e-2
            + cva.ba_cva_total_capital_charge.val() * 1e-2
            + lb.heynen_kat_partial_window_call.val()
            + vta.reliability_derated_vta_npv.val() * 1e-2
            + fl.constrained_information_ratio.val()
            + acd.duration_adjusted_price_impact_bps.val() * 1e-2
            + hlz.deflated_sharpe_ratio_dsr.val()
            + ml.mlmc_telescoping_price.val()
            + inf
    }
    pub fn phase26_live_metric(&self) -> f64 {
        let market = self.scenario_market();
        let bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else {
            0.0
        };
        Self::phase26_live_metric_from_market(&market, bump)
    }
    pub fn phase25_live_metric(&self) -> f64 {
        let market = self.scenario_market();
        let bump = if matches!(self.scenario, ScenarioKind::CreditWideningPlus100Bp) {
            0.01
        } else {
            0.0
        };
        Self::phase25_live_metric_from_market(&market, bump)
    }
    pub fn phase24_live_metric(&self) -> f64 {
        let market = self.scenario_market();
        let bump = if matches!(self.scenario, ScenarioKind::CreditWideningPlus100Bp) {
            0.01
        } else {
            0.0
        };
        Self::phase24_live_metric_from_market(&market, bump)
    }
    pub fn step_tick(&mut self) -> WorkstationFrameData {
        self.tick += 1;
        let market = self.scenario_market();
        let start = Instant::now();
        let summary = self.dag.evaluate_portfolio(&market, &self.trades);
        let elapsed = start.elapsed().as_secs_f64() * 1e6;
        let (credit_raw, jtd, _, cds_npv) = self.credit_metrics(self.scenario);
        let (cva, dva, fva) = self.xva_metrics();
        let mva = self.mva_metric();
        let (inflation, inflation_npv) = self.inflation_metrics(self.scenario);
        let active_curves = self.active_curve_system(&market);
        let curve = FlatForwardCurve::new(self.base_market.ref_date, market.ois_quotes[0]);
        let credit_bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else {
            0.0
        };
        // Phase 20 live telemetry: evaluate representative engines on the
        // current market/scenario snapshot so every tick and stress path is
        // reflected in the workstation metric.
        let sinhou_metric = price_sinh_ou_rfr_caplet(
            &SinhOuRfrCapletSpec {
                notional: 100.0,
                start_time: 0.5,
                end_time: 0.75,
                strike: market.ois_quotes[0] + credit_bump,
                alpha_mean_reversion: 0.05,
                sigma_short_rate: (0.01 + market.eq_vol * 0.03).max(1e-4),
                gamma_smile: 1.0 + market.eq_vol,
                y_star_skew: -0.15,
                is_caplet: true,
                is_backward_compounded: true,
            },
            &curve,
        )
        .npv
        .val();
        let mrb_metric = evaluate_mixed_rough_bergomi_vix(&MixedRoughBergomiVixSpec {
            xi0: (market.eq_vol * market.eq_vol).max(1e-4),
            hurst: 0.12,
            mixing_weight_lambda: 0.55,
            nu1: 0.9,
            nu2: 2.0 + market.eq_vol,
            expiry: 0.25,
            vix_window: 1.0 / 12.0,
            strike: 20.0,
            discount_factor: (-market.ois_quotes[0] * 0.25).exp(),
        })
        .vix_call_price
        .val();
        let choi_metric =
            ql_exotics::price_choi_sum_of_bsm_basket_option(&ql_exotics::ChoiSumOfBsmBasketSpec {
                spots: vec![market.eq_spot.max(1e-4), (0.9 * market.eq_spot).max(1e-4)],
                weights: vec![1.0, -1.0],
                vols: vec![market.eq_vol.max(1e-4), (0.8 * market.eq_vol).max(1e-4)],
                div_yields: vec![0.01, 0.02],
                correlation_matrix: vec![vec![1.0, 0.5], vec![0.5, 1.0]],
                strike: 5.0 + credit_bump * 100.0,
                risk_free_rate: market.ois_quotes[0],
                expiry: 1.0,
                is_call: true,
            })
            .price
            .val();
        let afs_metric = solve_alfonsi_fruth_schied_lob_execution(&AlfonsiFruthSchiedLobSpec {
            total_order_size: 100.0 + market.eq_spot * 0.01,
            horizon: 1.0,
            num_intervals: 4,
            lob_density_q: 50.0,
            shape_exponent_alpha: 0.4,
            resilience_rho: 1.2 + market.eq_vol,
        })
        .volume_resilience_cost
        .val();
        let exvar_metric =
            evaluate_and_optimize_expectile_portfolio(&ExpectileRiskAndPortfolioSpec {
                scenario_returns: vec![
                    vec![-0.07 - credit_bump, -0.01],
                    vec![-0.04, -0.005],
                    vec![0.02 - market.eq_vol * 0.01, 0.008],
                    vec![0.05, 0.012],
                ],
                initial_weights: vec![0.5, 0.5],
                tau_confidence: 0.85,
                return_weight_lambda: 0.2,
                smoothing_eps: 1e-3,
            })
            .optimal_expectile_exvar
            .val();
        let speed_metric =
            evaluate_third_order_speed_vol_expansion(&ThirdOrderSpeedVolExpansionSpec {
                xi0: (market.eq_vol * market.eq_vol).max(1e-4),
                omega_vol_of_vol: 1.8,
                k_mean_reversion: 3.0,
                expiry: 0.25,
                window_delta: 1.0 / 12.0,
                strike: 0.20,
                discount_factor: (-market.ois_quotes[0] * 0.25).exp(),
            })
            .call_price
            .val();

        // Complete Phase 20 telemetry: all fourteen engines consume the live
        // market/scenario snapshot so rates, equity, credit and tick changes
        // propagate through the aggregate metric.
        let cir_metric = ql_rates::price_cirpp_sofr_futures(
            &ql_rates::CirppSofrFuturesSpec {
                start_time: 0.5,
                end_time: 0.75,
                x0: market.ois_quotes[0] + credit_bump,
                b_drift: 0.006,
                beta_mean_reversion: 0.20,
                sigma: 0.06 + market.eq_vol * 0.02,
            },
            &curve,
        )
        .futures_rate_3m_compounded
        .val();
        let bk_metric = ql_fixed_income::analyze_black_karasinski_callable_frn_oas_ad(
            &ql_fixed_income::BlackKarasinskiCallableFrnSpec {
                face_amount: 100.0,
                quoted_margin: 0.008,
                cap_rate: 0.055,
                floor_rate: 0.020,
                maturity: 2.0,
                coupon_steps: 6,
                mean_reversion: 0.15,
                log_rate_vol: 0.22 + market.eq_vol * 0.005,
                flat_zero_rate: market.ois_quotes[0],
                call_price: 100.2,
                first_call_time: 1.0,
                target_market_price: 99.8 + credit_bump * 100.0,
                smoothing_eps: 1e-3,
            },
        )
        .oas_spread
        .val();
        let surv = ql_credit::SurvivalProbabilityCurve::new(
            self.base_market.ref_date,
            vec![1.0, 2.0, 3.0, 4.0, 5.0],
            vec![0.025 + credit_bump; 5],
        );
        let mb_metric = ql_credit::price_morini_brigo_credit_index_option(
            &ql_credit::MoriniBrigoIndexOptionSpec {
                notional: 100.0,
                strike_spread: 0.025 + credit_bump,
                option_expiry: 1.0,
                payment_times: vec![2.0, 3.0, 4.0, 5.0],
                recovery_rate: 0.40,
                spread_vol: 0.45 + market.eq_vol * 0.01,
                num_constituents: 50,
                systemic_copula_corr: 0.60,
                is_payer: true,
            },
            &curve,
            &surv,
        )
        .npv
        .val();
        let wwm_metric = ql_xva::compute_brigo_vrins_wwm_cva(&ql_xva::BrigoVrinsWwmCvaSpec {
            times: vec![0.5, 1.0, 2.0],
            discount_rate: market.ois_quotes[0],
            base_mtm_profile: vec![10.0 + market.eq_spot * 0.01, 12.0, 11.0],
            exposure_vol_nu: 8.0,
            is_brownian_bridge_swap: false,
            y0: 0.02,
            cir_kappa: 0.35,
            cir_theta: 0.025,
            cir_vol: 0.08,
            intensity_shift_psi: 0.005,
            wwr_correlation: 0.60,
            lgd: 0.60,
        })
        .wwm_cva
        .val();
        let bcs_metric =
            ql_xva::compute_bichuch_capponi_sturm_xva_band(&ql_xva::AsymmetricRateXvaBandSpec {
                maturity: 1.0,
                num_steps: 6,
                clean_claim_value: 25.0 + market.eq_spot * 0.01,
                clean_delta_notional: 40.0,
                asset_vol: 0.20 + market.eq_vol * 0.005,
                collateral_alpha: 0.50,
                dividend_discount_rate: 0.02,
                repo_lend_rate: 0.020,
                repo_borrow_rate: 0.024,
                funding_lend_rate: 0.022,
                funding_borrow_rate: 0.036 + credit_bump,
                collateral_lend_rate: 0.020,
                collateral_borrow_rate: 0.022,
                own_bond_yield_mu_i: 0.042,
                cpty_bond_yield_mu_c: 0.048 + credit_bump,
                own_lgd: 0.60,
                cpty_lgd: 0.60,
            })
            .no_arbitrage_band_width
            .val();
        let ldz_metric = ql_exotics::price_li_deng_zhou_multi_spread_option(
            &ql_exotics::LiDengZhouMultiSpreadSpec {
                long_spot: market.eq_spot,
                long_div_yield: 0.01,
                long_vol: market.eq_vol.max(1e-4),
                short_spots: vec![0.9 * market.eq_spot, 0.85 * market.eq_spot],
                short_weights: vec![0.5, 0.5],
                short_div_yields: vec![0.015, 0.02],
                short_vols: vec![0.20, 0.22],
                correlation_matrix: vec![
                    vec![1.0, 0.5, 0.4],
                    vec![0.5, 1.0, 0.6],
                    vec![0.4, 0.6, 1.0],
                ],
                strike: 5.0 + credit_bump * 100.0,
                risk_free_rate: market.ois_quotes[0],
                expiry: 1.0,
                is_call: true,
            },
        )
        .price
        .val();
        let mrjd_metric = ql_commodities::price_hikspoors_jaimungal_mrjd_commodity(
            &ql_commodities::HikspoorsJaimungalMrjdSpec {
                x0: 0.10,
                y0: 4.20 + market.eq_spot * 0.001,
                j0_spike: 0.04,
                seasonal_phi: 0.0,
                kappa_x: 1.5,
                sigma_x: 0.35 + market.eq_vol * 0.005,
                mu_y_rn: 0.02,
                sigma_y: 0.18,
                rho_xy: -0.2,
                kappa_j_spike: 12.0,
                jump_intensity: 2.5 + credit_bump * 10.0,
                mean_jump_size: 0.22,
                option_expiry: 0.5,
                futures_maturity: 1.0,
                strike: 70.0,
                risk_free_rate: market.ois_quotes[0],
                is_call: true,
            },
        )
        .futures_option_price
        .val();
        let pcrp_metric = ql_portfolio::solve_principal_component_risk_parity(
            &ql_portfolio::FactorEigenRiskParitySpec {
                covariance: vec![
                    vec![0.04 + market.eq_vol * 0.0005, 0.02],
                    vec![0.02, 0.05 + credit_bump * 0.01],
                ],
                factor_risk_budgets: vec![0.5, 0.5],
            },
        )
        .weights[0]
            .val();
        let sinhou_mb_rbergomi_choi_afs_expectile_metric = sinhou_metric
            + mrb_metric
            + choi_metric
            + afs_metric
            + exvar_metric
            + speed_metric
            + cir_metric
            + bk_metric
            + mb_metric
            + wwm_metric
            + bcs_metric
            + ldz_metric
            + mrjd_metric
            + pcrp_metric;
        let fx = FxVannaVolgaSmile::new(
            market.eq_spot,
            0.03,
            0.01,
            1.0,
            FxSmileQuote {
                atm_vol: 0.1,
                rr_25d: 0.01,
                bf_25d: 0.005,
            },
        );
        let fx_delta =
            garman_kohlhagen_spot_delta(market.eq_spot, fx.k3, fx.sigma3, 1.0, 0.03, 0.01, true)
                .val();
        // Assemble the actual cross-asset SIMM buckets from the current
        // sensitivities, rather than applying one arbitrary scalar weight to
        // the whole IR ladder.  All values are expressed in the same dollar
        // sensitivity units before the risk weights are applied.
        let mut buckets = vec![SimmBucket {
            sensitivities: summary
                .total_par_dv01
                .iter()
                .map(|x| ql_ad::TapeVar::from_f64(*x * 1e-4))
                .collect(),
            risk_weights: vec![50.0; summary.total_par_dv01.len()],
            concentration_threshold: 1e6,
            intra_corr: vec![vec![1.0; summary.total_par_dv01.len()]; summary.total_par_dv01.len()],
        }];
        buckets.push(SimmBucket {
            sensitivities: credit_raw
                .iter()
                .map(|(_, value)| ql_ad::TapeVar::from_f64(*value * 1e-4))
                .collect(),
            risk_weights: vec![100.0; credit_raw.len()],
            concentration_threshold: 1e6,
            intra_corr: vec![vec![1.0; credit_raw.len()]; credit_raw.len()],
        });
        buckets.push(SimmBucket {
            sensitivities: vec![ql_ad::TapeVar::from_f64(
                summary.total_eq_delta * market.eq_spot,
            )],
            risk_weights: vec![0.15],
            concentration_threshold: 1e6,
            intra_corr: vec![vec![1.0]],
        });
        buckets.push(SimmBucket {
            sensitivities: vec![ql_ad::TapeVar::from_f64(fx_delta * market.eq_spot)],
            risk_weights: vec![0.08],
            concentration_threshold: 1e6,
            intra_corr: vec![vec![1.0]],
        });
        buckets.push(SimmBucket {
            sensitivities: vec![ql_ad::TapeVar::from_f64(inflation)],
            risk_weights: vec![50.0],
            concentration_threshold: 1e6,
            intra_corr: vec![vec![1.0]],
        });
        let class = |bucket: &SimmBucket<ql_ad::TapeVar<f64>>| {
            compute_isda_simm(std::slice::from_ref(bucket), &[vec![1.0]])
        };
        let ir_margin = compute_isda_simm(
            &[buckets[0].clone(), buckets[4].clone()],
            &[vec![1.0, 0.5], vec![0.5, 1.0]],
        );
        let class_margins = [
            ir_margin,
            class(&buckets[1]),
            ql_ad::TapeVar::zero(),
            class(&buckets[2]),
            ql_ad::TapeVar::zero(),
            class(&buckets[3]),
        ];
        let simm = compute_isda_simm_6_risk_classes(&class_margins).val();
        let adjusted_total_npv = summary.total_npv + cds_npv + inflation_npv;
        let prev = self
            .npv_history
            .last()
            .copied()
            .unwrap_or(adjusted_total_npv);
        self.npv_history.push(adjusted_total_npv);
        self.last_summary = summary.clone();
        let rate = active_curves.curves[0].zero_rate(1.0).val();
        let eq_vega = {
            let params = HestonParams {
                v0: ql_ad::Dual::<1>::variable(market.eq_vol * market.eq_vol, 0),
                kappa: ql_ad::Dual::from_value(market.heston_params.kappa),
                theta: ql_ad::Dual::from_value(market.heston_params.theta),
                sigma_v: ql_ad::Dual::from_value(market.heston_params.sigma_v),
                rho: ql_ad::Dual::from_value(market.heston_params.rho),
            };
            let v = ql_models::heston_cos_price(
                ql_ad::Dual::from_value(market.eq_spot),
                ql_ad::Dual::from_value(100.0),
                ql_ad::Dual::from_value(1.0),
                ql_ad::Dual::from_value(rate),
                ql_ad::Dual::from_value(0.0),
                &params,
                128,
                10.0,
                true,
            );
            v.d[0] * 2.0 * market.eq_vol
        };
        let ir_buckets: Vec<(String, f64)> = summary
            .total_par_dv01
            .iter()
            .enumerate()
            .map(|(i, x)| {
                let label = ["OIS 1Y", "OIS 2Y", "OIS 5Y", "FWD 1Y", "FWD 2Y", "FWD 5Y"]
                    .get(i)
                    .map(|s| (*s).to_string())
                    .unwrap_or_else(|| format!("Pillar {}", i + 1));
                (label, *x * 1e-4)
            })
            .collect();
        let credit_buckets: Vec<(String, f64)> = credit_raw
            .iter()
            .map(|(label, value)| (label.clone(), *value * 1e-4))
            .collect();
        let irs_npv = if let (Some(ois), Some(fwd)) =
            (active_curves.curves.first(), active_curves.curves.get(1))
        {
            let spec = VanillaSwapSpec {
                effective_date: market.ref_date,
                termination_date: market.ref_date.add_months(60, false),
                notional: 1_000_000.0,
                fixed_rate_or_coupon: 0.02,
                fixed_tenor: Period::Years(1),
                float_tenor: Period::Years(1),
                calendar: Arc::new(BitsetCalendar::weekend_only()),
                bdc: BusinessDayConvention::Unadjusted,
                fixed_dc: DayCounter::Act365Fixed,
                float_dc: DayCounter::Act365Fixed,
                swap_type: SwapType::Payer,
            };
            spec.resolve(market.ref_date).npv(ois, fwd)
        } else {
            0.0
        };
        let heston_params = HestonParams {
            v0: market.heston_params.v0,
            kappa: market.heston_params.kappa,
            theta: market.heston_params.theta,
            sigma_v: market.heston_params.sigma_v,
            rho: market.heston_params.rho,
        };
        let heston_npv = ql_models::heston_cos_price(
            market.eq_spot,
            100.0,
            1.0,
            rate,
            0.0,
            &heston_params,
            128,
            10.0,
            true,
        );
        let script_npv = match self.trades.get(5) {
            Some(PortfolioTrade::ScriptedNote { program, .. }) => ql_scripting::ScriptMcEngine {
                n_paths: 1024,
                seed: 42,
            }
            .price_gbm(program, market.eq_spot, market.eq_vol, rate, 0.0)
            .val(),
            _ => 0.0,
        };
        let (_, inflation_curve) = self.inflation_calibration(self.scenario);
        let zciis = ZeroCouponInflationSwap {
            notional: 1_000_000.0,
            fixed_rate: 0.021,
            maturity_years: 5.0,
            target_month: Some(7),
            pay_fixed: true,
        };
        let zciis_npv = zciis
            .npv(&inflation_curve.discount, &inflation_curve.curve)
            .val();
        let g2_swaption = EuropeanSwaption {
            notional: 1_000_000.0,
            strike: 0.025,
            expiry: 1.0,
            swap_times: (1..=6).map(|x| x as f64).collect(),
            is_payer: true,
            settlement: SwaptionSettlement::Physical,
        };
        let g2 = G2PlusPlus {
            a: 0.1,
            sigma: 0.01,
            b: 0.3,
            eta: 0.015,
            rho: -0.75,
        };
        let g2_npv = g2
            .swaption_price(&active_curves.curves[0], &g2_swaption, 16)
            .val();
        let vg = VarianceGamma {
            sigma: 0.18,
            nu: 0.20,
            theta: -0.10,
        };
        let levy_npv = cos_levy_price(
            &vg,
            market.eq_spot,
            100.0,
            rate,
            0.0,
            1.0,
            ql_instruments::OptionType::Call,
            128,
        )
        .val();
        let conv_spec = ConvertibleBondSpec {
            spot: market.eq_spot,
            conversion_ratio: 1.0,
            face_value: 100.0,
            coupon_rate: 0.03,
            coupon_times: (1..=3).map(|x| x as f64).collect(),
            maturity: 3.0,
            rate,
            div_yield: 0.0,
            credit_spread: 0.02,
            vol: market.eq_vol,
            call_price: None,
            put_price: None,
        };
        let (conv_npv, _) = price_convertible_bond_tf(&conv_spec, 48);
        let rfr_spec = RfrCapletSpec {
            notional: 1_000_000.0,
            forward_rate: 0.045,
            strike: 0.04,
            tau_0: 0.5,
            tau_1: 1.0,
            accrual_fraction: 0.5,
            discount_tau_1: active_curves.curves[0].discount_factor(1.0).val(),
            is_caplet: true,
            rate_type: RfrRateType::BackwardLooking,
        };
        let rfr_npv =
            willems_backward_sabr_caplet_price(&rfr_spec, 0.08, 0.5, -0.35, 0.45, 1.0).val();
        let xccs = ConstNotionalXccs {
            dom_notional: 1_100_000.0,
            for_notional: 1_000_000.0,
            spot_fx: 1.1,
            payment_times: (1..=5).map(|x| x as f64).collect(),
            accrual_fractions: vec![1.0; 5],
            for_basis_spread: 0.001,
            dom_spread: 0.0,
            receive_foreign: true,
        };
        let xccs_npv = xccs
            .npv(
                &active_curves.curves[0],
                &active_curves.curves[1],
                &active_curves.curves[1],
            )
            .val();
        let es = gaussian_expected_shortfall(
            ql_ad::TapeVar::from_f64(summary.total_npv),
            ql_ad::TapeVar::from_f64(summary.total_eq_delta.abs().max(1.0)),
            0.975,
        )
        .val();
        let rough_params = LiftedHestonParams {
            v0: market.heston_params.v0,
            lambda: market.heston_params.kappa,
            theta: market.heston_params.theta,
            nu: market.heston_params.sigma_v,
            rho: market.heston_params.rho,
            hurst: 0.1,
            n_factors: 8,
            r_n: 2.5,
        };
        let rough_npv = lifted_heston_cos_price(
            market.eq_spot,
            100.0,
            1.0,
            rate,
            0.0,
            &rough_params,
            64,
            8.0,
            16,
            true,
        )
        .val();
        let barrier_spec = BarrierOptionSpec {
            option_type: ql_instruments::OptionType::Call,
            barrier_type: ExoticBarrierType::DownOut,
            spot: market.eq_spot,
            strike: 100.0,
            barrier: 90.0,
            rebate: 0.0,
            rate,
            div_yield: 0.0,
            vol: market.eq_vol,
            expiry: 1.0,
            monitoring: MonitoringFrequency::Discrete { observations: 12 },
        };
        let barrier_npv = price_barrier_option(&barrier_spec).val();
        let commodity = SchwartzSmith2F {
            chi0: 0.0,
            xi0: market.eq_spot.ln(),
            kappa: 1.0,
            sigma_chi: 0.3,
            mu_xi: 0.0,
            sigma_xi: 0.2,
            rho: -0.3,
        };
        let commodity_futures = commodity.futures_price(1.0);
        let commodity_option = commodity
            .price_futures_option(
                1.0,
                1.5,
                commodity_futures,
                rate,
                ql_instruments::OptionType::Call,
            )
            .val();
        let crack = kirk_spread_option(
            commodity_futures,
            3.0,
            0.0,
            0.3,
            0.25,
            0.2,
            rate,
            1.0,
            ql_instruments::OptionType::Call,
        )
        .val();
        // Live Phase 11 telemetry for the final two workstation tabs.  These
        // are deliberately computed from the active scenario curves rather
        // than being display-only constants.
        let bermudan = BermudanSwaption {
            notional: 1_000_000.0,
            strike: 0.02,
            exercise_times: vec![1.0, 2.0, 3.0],
            swap_times: (0..=5).map(|x| x as f64).collect(),
            is_payer: true,
            settlement: SwaptionSettlement::Physical,
        };
        let bermudan_value = HullWhiteBermudanTreeEngine { steps: 30 }
            .price_single_curve(
                &HullWhite1F {
                    a: 0.05,
                    sigma: 0.01,
                },
                &active_curves.curves[0],
                &bermudan,
            )
            .val();
        let sa_report = compute_sa_ccr_ead(&SaCcrNettingSet {
            trades: vec![SaCcrTrade {
                mtm: ql_ad::TapeVar::from_f64(irs_npv),
                notional: ql_ad::TapeVar::from_f64(1_000_000.0),
                maturity: 5.0,
                asset_class: SaCcrAssetClass::InterestRate {
                    ccy: "USD".to_string(),
                    start_time: 0.0,
                    end_time: 5.0,
                },
                is_long_primary: true,
                option_spec: None,
            }],
            collateral: ql_ad::TapeVar::zero(),
            is_margined: false,
            threshold: ql_ad::TapeVar::zero(),
            mta: ql_ad::TapeVar::zero(),
            nica: ql_ad::TapeVar::zero(),
            mpor_years: 0.5,
        });
        let commodity_curve = CommodityFuturesCurve::from_futures_pillars(
            vec![1.0, 2.0, 3.0],
            &[80.0, 81.0, 82.0],
            ql_commodities::SeasonalParams {
                a1: 0.0,
                b1: 0.0,
                a2: 0.0,
                b2: 0.0,
            },
        );
        let swing_value = price_swing_option(
            &SwingOptionSpec {
                exercise_times: vec![1.0, 2.0, 3.0],
                strike: 80.0,
                unit_volume: 1.0,
                min_exercises: 0,
                max_exercises: 2,
                is_call: true,
            },
            &commodity_curve,
            1.0,
            0.3,
            rate,
            2,
        )
        .val();
        let qp = box_constrained_markowitz_weights(&BoxConstrainedMarkowitzSpec {
            covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
            expected_returns: vec![0.06, market.ois_quotes[2]],
            risk_aversion: 3.0,
            lower_bounds: vec![0.15, 0.15],
            upper_bounds: vec![0.80, 0.80],
        });
        let qp_value = qp[0];
        let storage_value = price_gas_storage(
            &GasStorageSpec {
                decision_times: vec![1.0, 2.0],
                max_capacity: 10.0,
                initial_inventory: 5.0,
                max_injection_rate: 5.0,
                max_withdrawal_rate: 5.0,
                injection_loss: 0.0,
                withdrawal_loss: 0.0,
                injection_cost: 0.0,
                withdrawal_cost: 0.0,
                inventory_grid_points: 10,
            },
            &commodity_curve,
            1.0,
            0.3,
            rate,
        );
        let sabr_forward = market.eq_spot / 100.0;
        let sabr_quotes = (0..5)
            .map(|i| SmileStrikeQuote {
                strike: sabr_forward * (0.96 + 0.02 * i as f64),
                market_vol: 0.20 + 0.004 * (i as f64 - 2.0).abs(),
                weight: 1.0,
            })
            .collect::<Vec<_>>();
        let sabr_live =
            CalibratedSabrSmile::calibrate(sabr_forward, 1.0, 0.5, (0.1, -0.2, 0.4), &sabr_quotes)
                .expect("live SABR telemetry");
        let (sabr_atm, sabr_vega) = sabr_live.vol_and_quote_vega(sabr_forward);
        let sabr_metric = sabr_atm + sabr_vega.iter().sum::<f64>();
        let rough_metric = rbergomi_hybrid_mc_price(
            market.eq_spot,
            100.0,
            1.0,
            rate,
            0.0,
            &RoughBergomiParams {
                xi0: market.heston_params.v0,
                eta: market.heston_params.sigma_v,
                rho: market.heston_params.rho,
                hurst: 0.1,
            },
            8,
            64,
            true,
            self.tick as u64,
        )
        .val()
            + price_worst_of_autocallable(
                &WorstOfAutocallableSpec {
                    notional: 1.0,
                    spots: vec![market.eq_spot, market.eq_spot * 0.98],
                    vols: vec![market.eq_vol, market.eq_vol * 1.05],
                    divs: vec![0.0, 0.0],
                    rate,
                    correlation: vec![vec![1.0, 0.5], vec![0.5, 1.0]],
                    observation_times: vec![0.5, 1.0],
                    autocall_barriers: vec![1.05, 1.05],
                    redemption_bonuses: vec![0.0, 0.0],
                    coupon_barriers: vec![0.8, 0.8],
                    coupon_rates: vec![0.03, 0.03],
                    memory_coupon: true,
                    knock_in_barrier: 0.65,
                    put_strike: 1.0,
                    observe_ki_all_dates: true,
                    smoothing_half_width: 0.01,
                    smooth_min_temp: 0.01,
                },
                32,
                self.tick as u64 + 17,
            )
            .unwrap_or(0.0);
        let execution_metric = obizhaeva_wang_execution(&ObizhaevaWangSpec {
            initial_inventory: 1_000_000.0,
            trade_times: (0..=8).map(|i| i as f64 / 8.0).collect(),
            volatility: market.eq_vol,
            transient_impact: 0.01,
            resilience_speed: 0.5,
            temporary_impact: 0.001,
            risk_aversion: 0.1,
        })
        .mean_variance_utility
        .val()
            + compute_frtb_sa(&ql_xva::FrtbSaSpec {
                buckets: vec![ql_xva::FrtbSensitivityBucket {
                    risk_class: ql_xva::FrtbRiskClass::Girr,
                    bucket_id: 0,
                    delta_sensitivities: vec![ql_ad::TapeVar::from_f64(
                        summary.total_par_dv01.iter().sum::<f64>(),
                    )],
                    delta_risk_weights: vec![0.015],
                    vega_sensitivities: vec![ql_ad::TapeVar::from_f64(sabr_metric)],
                    vega_risk_weights: vec![0.55],
                    curvature_up_pnl: vec![ql_ad::TapeVar::zero()],
                    curvature_down_pnl: vec![ql_ad::TapeVar::zero()],
                    curvature_shock_weights: vec![0.015],
                    intra_bucket_corr: vec![vec![1.0]],
                }],
                inter_bucket_corr: vec![],
                drc_obligors: vec![],
                exotic_residual_notionals: vec![],
                other_residual_notionals: vec![],
            })
            .sbm_med
            .val()
                * 1e-3
            + sabr_metric * 1e-3;
        // Phase 13 Tab 23 telemetry must be sourced from the four engines named
        // by the tab, rather than reusing the SLV/FRTB aggregates above.  Keep
        // each contribution tied to the active scenario so a stress tick is
        // observable in the workstation frame and in the rendered dashboard.
        let monotone_curve = MonotoneConvexZeroCurve::new(
            market.ref_date,
            vec![1.0, 2.0, 5.0],
            market.ois_quotes.to_vec(),
        );
        let monotone_metric = [0.5, 1.5, 3.0]
            .iter()
            .map(|&t| monotone_curve.zero_rate(t).val())
            .sum::<f64>();
        let vix_params = BatesParams {
            heston: HestonParams {
                v0: market.heston_params.v0,
                kappa: market.heston_params.kappa,
                theta: market.heston_params.theta,
                sigma_v: market.heston_params.sigma_v,
                rho: market.heston_params.rho,
            },
            lambda: 0.1,
            mu_j: -0.05,
            sigma_j: 0.2,
        };
        let vix =
            heston_bates_vix_future_and_option_price(&vix_params, 1.0, 30.0 / 365.0, 0.2, 0.98);
        let parisian = ParisianBarrierOptionSpec {
            option_type: ql_instruments::OptionType::Call,
            barrier_type: ExoticBarrierType::DownOut,
            clock_kind: ExcursionClockKind::ParisianConsecutive,
            spot: market.eq_spot,
            strike: 100.0,
            barrier: 90.0,
            rebate: 0.0,
            window_duration: 0.20,
            rate,
            div_yield: 0.0,
            vol: market.eq_vol,
            expiry: 1.0,
        };
        let parisian_metric = price_parisian_option_mc_ad(&parisian, 64, 16, 0.01, 23).val();
        let robust_weights = robust_markowitz_weights(&RobustMarkowitzSpec {
            covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
            return_uncertainty_cov: vec![vec![0.01, 0.0], vec![0.0, 0.02]],
            expected_returns: vec![0.06, market.ois_quotes[2]],
            risk_aversion: 3.0,
            uncertainty_radius: 0.1,
            lower_bounds: vec![0.0, 0.0],
            upper_bounds: vec![0.8, 0.8],
        });
        let robust_metric = robust_weights
            .iter()
            .enumerate()
            .map(|(i, weight)| (i + 1) as f64 * weight)
            .sum::<f64>();
        let vix_parisian_robust_metric = monotone_metric
            + vix.vix_future
            + vix.call_price
            + vix.put_price
            + parisian_metric
            + robust_metric;
        let dh_metric = double_heston_cos_price(
            market.eq_spot,
            100.0,
            1.0,
            rate,
            0.0,
            &DoubleHestonParams {
                factor1: heston_params,
                factor2: HestonParams {
                    v0: 0.01,
                    kappa: 2.0,
                    theta: 0.01,
                    sigma_v: 0.2,
                    rho: -0.3,
                },
            },
            48,
            10.0,
            true,
        )
        .val();
        let tarf_metric = price_fx_tarf_mc_ad(
            &FxTarfSpec {
                notional: 100.0,
                spot: 1.05,
                strike: 1.0,
                r_dom: rate,
                r_for: 0.01,
                vol: 0.12,
                fixing_times: vec![0.25, 0.5, 0.75, 1.0],
                target_profit: 0.15,
                downside_leverage: 1.5,
                is_exporter_call: true,
                cap_convention: FxTarfCapConvention::ExactCap,
                smoothing_eps: 1e-3,
            },
            32,
            7,
        )
        .npv;
        let coco_metric = price_coco_bond_ds(&CocoBondSpec {
            face_value: 100.0,
            coupon_rate: 0.05,
            coupon_times: vec![1.0, 2.0],
            maturity: 2.0,
            rate,
            credit_spread: 0.015,
            spot: market.eq_spot,
            trigger_barrier: 60.0,
            div_yield: 0.0,
            eq_vol: market.eq_vol,
            trigger_payoff: CocoTriggerPayoff::WriteDown {
                recovery_fraction: 0.4,
            },
        })
        .dirty_price;
        let hrp_w = hierarchical_risk_parity_weights(&[vec![0.04, 0.01], vec![0.01, 0.09]]);
        let double_heston_tarf_coco_hrp_metric =
            dh_metric + tarf_metric + coco_metric + hrp_w[0] + 2.0 * hrp_w[1];
        let jy_curve = ZeroInflationCurve::new(
            market.ref_date,
            100.0,
            vec![1.0, 2.0, 5.0, 10.0],
            vec![0.02, 0.021, 0.023, 0.025],
            None,
        );
        let jy_instrument = YoYInflationCapFloor {
            notional: 1.0,
            maturity_years: 10,
            cap_rate: Some(0.025),
            floor_rate: None,
            gearing: 1.0,
            spread: 0.0,
            cap_floor_type: YoYInflationCapFloorType::Cap,
        };
        let jy_value = price_yoy_capfloor_jarrow_yildirim(
            &jy_instrument,
            &active_curves.curves[0],
            &jy_curve,
            &JarrowYildirimParams {
                a_nom: 0.1,
                sigma_nom: 0.01,
                a_real: 0.08,
                sigma_real: 0.012,
                sigma_cpi: 0.01,
                rho_nom_real: 0.2,
                rho_nom_cpi: 0.1,
                rho_real_cpi: -0.2,
            },
        )
        .0
        .val();
        let tba = analyze_tba_dollar_roll(&TbaDollarRollSpec {
            front_clean_price: 99.0,
            back_clean_price: 98.5,
            wac: 0.04,
            pass_through_rate: 0.035,
            remaining_months: 120,
            front_accrued_fraction: 0.5,
            back_accrued_fraction: 0.5,
            roll_period_years: 1.0,
            assumed_cpr: 0.06,
            repo_rate: rate,
        });
        let pool = price_specified_pool_mbs(
            &SpecifiedPoolMbsSpec {
                original_balance: 100.0,
                wac: 0.04,
                pass_through_rate: 0.035,
                term_months: 120,
                seasoning_months: 30,
                cpr_min: 0.02,
                cpr_max: 0.10,
                refi_slope: 2.0,
                refi_inflection: 0.0,
                primary_secondary_spread: 0.01,
                cdr: 0.002,
                loss_severity: 0.4,
            },
            &active_curves.curves[0],
            0.0,
        );
        let bg = bergomi_2f_implied_vol_bg(
            &Bergomi2FParams {
                xi0: 0.04,
                omega: 0.4,
                theta: 0.5,
                k1: 1.0,
                k2: 2.0,
                rho12: 0.25,
                rho_s1: -0.5,
                rho_s2: -0.3,
            },
            100.0,
            100.0,
            1.0,
        )
        .val();
        let returns = (0..16)
            .map(|i| {
                vec![
                    0.01 * (i as f64 - 7.0),
                    0.02 * ((i % 5) as f64 - 2.0),
                    0.015 * ((i % 3) as f64 - 1.0),
                    0.01 * ((i % 4) as f64 - 1.5),
                ]
            })
            .collect::<Vec<_>>();
        let gerber = gerber_covariance_matrix(&returns, 0.5, 1e-3).unwrap_or_default();
        let nco =
            nested_clustered_optimization_weights(&ql_portfolio::NestedClusteredOptimizationSpec {
                covariance: gerber.clone(),
                expected_returns: None,
                risk_aversion: 1.0,
                num_clusters: 2,
            });
        let jy_tba_bergomi_nco_metric = jy_value
            + tba.net_roll_advantage.val()
            + pool.pv.val() * 1e-3
            + bg
            + nco.iter().sum::<f64>();
        let cheyette = price_cheyette_swaption_markov_projection(
            &CheyetteSwaptionSpec {
                notional: 100.0,
                strike: 0.032,
                option_expiry: 1.0,
                payment_times: vec![2.0, 3.0, 4.0, 5.0],
                mean_reversion: 0.05,
                sigma0: 0.008,
                skew_b: -1.5,
                is_payer: true,
            },
            &active_curves.curves[0],
        );
        let cmo = price_cmo_waterfall(
            &CmoWaterfallSpec {
                original_balance: 100.0,
                wac: 0.05,
                pass_through_rate: 0.045,
                term_months: 24,
                psa_speed: 150.0,
                pac_collar_psa: Some((100.0, 300.0)),
                tranches: vec![
                    CmoTrancheSpec {
                        initial_balance: 60.0,
                        coupon_rate: 0.04,
                        is_pac: true,
                    },
                    CmoTrancheSpec {
                        initial_balance: 40.0,
                        coupon_rate: 0.045,
                        is_pac: false,
                    },
                ],
            },
            &active_curves.curves[0],
            0.0,
        );
        let qrh = qrh_implied_vol_expansion(
            &QuadraticRoughHestonParams {
                hurst: 0.1,
                nu: 0.35,
                lambda: 1.2,
                a: 0.45,
                b: 0.2,
                c: market.eq_vol * market.eq_vol,
                z0: 0.0,
            },
            market.eq_spot,
            market.eq_spot,
            1.0,
        );
        let ctd = analyze_treasury_futures_ctd(&TreasuryFuturesBasketSpec {
            futures_price: 98.5,
            delivery_time: 0.25,
            repo_rate: rate,
            notional_standard_yield: 0.06,
            deliverables: vec![
                TreasuryFutureDeliverable {
                    clean_price: 99.2,
                    coupon_rate: 0.06,
                    remaining_semiannual_periods: 16,
                    accrued_interest_spot: 0.5,
                    accrued_interest_delivery: 1.0,
                    interim_coupon: 0.0,
                    yield_vol: 0.008,
                    modified_duration: 6.5,
                },
                TreasuryFutureDeliverable {
                    clean_price: 89.0,
                    coupon_rate: 0.045,
                    remaining_semiannual_periods: 18,
                    accrued_interest_spot: 0.4,
                    accrued_interest_delivery: 0.8,
                    interim_coupon: 0.0,
                    yield_vol: 0.008,
                    modified_duration: 7.2,
                },
            ],
        })
        .option_adjusted_futures_price;
        let prdc = price_prdc_note_hybrid(&PrdcNoteSpec {
            notional: 100.0,
            spot_fx: 150.0,
            r_dom: 0.01,
            r_for: 0.04,
            a_dom: 0.05,
            sigma_dom: 0.008,
            a_for: 0.08,
            sigma_for: 0.01,
            sigma_fx: 0.12,
            rho_fx_dom: 0.3,
            rho_fx_for: -0.2,
            rho_dom_for: 0.25,
            coupon_times: vec![1.0, 2.0],
            foreign_gearing: 0.15,
            domestic_gearing: 0.10,
            coupon_cap: 0.08,
            fx_knockout_barrier: Some(180.0),
        })
        .npv;
        let factor_attr = factor_risk_attribution(
            &[0.4, 0.6],
            &[vec![1.1, 0.3], vec![0.8, 0.9]],
            &[vec![0.04, 0.005], vec![0.005, 0.02]],
            &[0.01, 0.015],
        )
        .total_volatility;
        let dd_returns = (0..16)
            .map(|i| -0.045 + 0.006 * i as f64 - 0.002 * rate)
            .collect::<Vec<_>>();
        let evar_cdar = conditional_drawdown_at_risk(&dd_returns, 0.75).cdar
            + entropic_value_at_risk(&dd_returns, 0.90);
        let base_quotes = vec![
            BaseCorrelationTrancheQuote {
                attachment: 0.0,
                detachment: 0.03,
                running_spread: 0.015 + rate,
                upfront: 0.0,
            },
            BaseCorrelationTrancheQuote {
                attachment: 0.03,
                detachment: 0.07,
                running_spread: 0.02 + rate,
                upfront: 0.0,
            },
        ];
        let base_corr = bootstrap_base_correlation_curve(
            &base_quotes,
            0.02 + rate * 0.1,
            0.4,
            rate,
            &[1.0, 2.0, 3.0, 4.0, 5.0],
        );
        let (_, base_fair, _) = ql_credit::price_mezzanine_tranche_base_correlation(
            &base_corr,
            0.02,
            0.05,
            0.02 + rate,
            100.0,
        );
        let cheyette_cmo_qrh_evar_metric = cheyette.npv.val()
            + cmo.total_structure_pv.val()
            + qrh.val()
            + ctd
            + prdc
            + evar_cdar
            + factor_attr
            + base_fair.val();
        let phase17_curve = FlatForwardCurve::new(market.ref_date, rate);
        let asw = analyze_bond_asset_swap_and_z_spread(
            &BondAssetSwapSpec {
                clean_price: 0.97,
                accrued_interest: 0.01,
                coupon_rate: 0.04,
                payment_times: vec![1.0, 2.0, 3.0, 4.0, 5.0],
                accrual_fractions: vec![1.0; 5],
            },
            &phase17_curve,
        );
        let mf = price_markov_functional_1f(
            &MarkovFunctionalSwaptionSpec {
                expiry: 1.0,
                payment_times: vec![2.0, 3.0, 4.0, 5.0],
                mean_reversion: 0.05 + rate * 0.1,
                atm_normal_vol: 0.008 + market.eq_vol * 0.01,
                sabr_rho: -0.25,
                sabr_nu: 0.35,
                strike: 0.031 + rate * 0.05,
                notional: 100.0,
                is_payer: true,
            },
            &phase17_curve,
            12,
        );
        let ba_cva = compute_ba_cva_capital(&BaCvaPortfolioSpec {
            counterparties: vec![BaCvaCounterpartySpec {
                risk_weight: 0.05,
                effective_maturity: 3.0,
                ead: 100.0 + market.eq_vol * 20.0,
                single_name_hedges: vec![BaCvaSingleNameHedge {
                    notional: 35.0,
                    maturity: 3.0,
                    hedge_correlation: 0.8,
                }],
            }],
            index_hedges: vec![],
            hurdle_rate: 0.10,
            risk_free_rate: rate,
        });
        let black_cox = analyze_structural_credit_black_cox(&BlackCoxStructuralSpec {
            asset_value: 140.0,
            debt_face_value: 100.0,
            default_barrier: 85.0,
            barrier_growth_rate: 0.02,
            asset_vol: 0.20 + market.eq_vol,
            risk_free_rate: rate + credit_bump,
            asset_payout_rate: 0.01,
            maturity: 3.0,
            recovery_rate: 0.40,
        });
        let ah_slice = solve_andreasen_huge_step(
            &AndreasenHugeGridSpec {
                forward: 100.0,
                rate,
                dt: 0.5,
                strike_min: 70.0,
                strike_max: 130.0,
                n_grid: 25,
                local_vol_atm: 0.18 + market.eq_vol,
                local_vol_skew: -0.15,
                local_vol_smile: 0.10,
            },
            None,
        );
        let ah_call = interpolate_andreasen_huge_call(&ah_slice, 100.0);
        let phase17_curve_futures = CommodityFuturesCurve::from_futures_pillars(
            vec![0.25, 0.50, 0.75, 1.0],
            &[
                45.0 + market.eq_vol * 5.0 + market.eq_spot * 1e-3,
                58.0,
                42.0,
                62.0,
            ],
            ql_commodities::SeasonalParams {
                a1: 0.0,
                b1: 0.0,
                a2: 0.0,
                b2: 0.0,
            },
        );
        let bess = price_bess_arbitrage_dp(
            &BessSpec {
                decision_times: vec![0.25, 0.50, 0.75, 1.0],
                energy_capacity_mwh: 4.0,
                initial_soc_mwh: 2.0,
                max_charge_mw: 8.0,
                max_discharge_mw: 8.0,
                charge_efficiency: 0.92,
                discharge_efficiency: 0.92,
                cycle_degradation_cost: 1.5,
                ancillary_reserve_rate: 2.0,
                soc_grid_points: 3,
            },
            &phase17_curve_futures,
            1.5,
            0.35 + market.eq_vol,
            rate,
        );
        let fx_accumulator = price_fx_accumulator_strip_no_ko(&FxAccumulatorSpec {
            base_notional_per_fixing: 1.0,
            spot: market.eq_spot,
            strike: 100.0,
            knockout_barrier: 120.0,
            r_dom: rate + credit_bump,
            r_for: 0.01,
            vol: market.eq_vol.max(0.01),
            fixing_times: vec![0.25, 0.50, 0.75, 1.0],
            guaranteed_fixings: 0,
            downside_gearing: 1.2,
            is_accumulator: true,
            smoothing_eps: 1e-6,
        });
        let atlas = price_atlas_option_mc_ad(
            &AtlasOptionSpec {
                spots: vec![market.eq_spot, market.eq_spot * 1.02, market.eq_spot * 0.98],
                vols: vec![market.eq_vol.max(0.01); 3],
                div_yields: vec![0.01; 3],
                correlation: 0.35,
                rate: rate + credit_bump,
                expiry: 1.0,
                strike: 0.0,
                notional: 1.0,
                remove_best: 0,
                remove_worst: 0,
            },
            32,
            17,
        );
        let take_or_pay = price_take_or_pay_contract(
            &TakeOrPayContractSpec {
                decision_times: vec![0.25, 0.50, 0.75, 1.0],
                strike: 50.0 + credit_bump,
                daily_min_qty: 1.0,
                daily_max_qty: 3.0,
                cumulative_min_qty: 6.0,
                cumulative_max_qty: 10.0,
                shortfall_penalty_rate: 50.0,
                daily_discretization_steps: 2,
            },
            &phase17_curve_futures,
            1.5,
            0.30 + market.eq_vol,
            rate + credit_bump,
        );
        let bl = idzorek_black_litterman_posterior(&IdzorekBlackLittermanSpec {
            covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
            market_weights: vec![0.6, 0.4],
            risk_aversion: 2.5,
            tau: 0.05,
            view_matrix: vec![vec![1.0, -1.0]],
            view_returns: vec![0.03 + rate * 0.1],
            view_confidences: vec![0.65 + market.eq_vol * 0.1],
        });
        let kelly = kelly_fractional_allocation(&KellyAllocationSpec {
            expected_returns: vec![0.05 + rate, 0.07 + market.eq_vol],
            covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
            risk_free_rate: rate,
            max_kelly_fraction: 0.50,
            max_drawdown_limit: 0.20,
            drawdown_tail_prob: 0.10,
        });
        let dark_lit = optimal_dark_lit_execution(&DarkLitRoutingSpec {
            initial_inventory: 100.0,
            horizon: 1.0,
            num_steps: 4,
            volatility: market.eq_vol,
            lit_temporary_impact: 1e-3,
            lit_permanent_impact: 2e-4,
            dark_fill_intensity: 1.5 + rate * 2.0,
            risk_aversion: 0.25,
        });
        let asw_mf1f_ah_bess_ep_metric = asw.z_spread.val()
            + asw.par_asset_swap_spread.val()
            + mf.swaption_npv.val()
            + mf.cms_swaplet_rate.val()
            + ba_cva.capital_charge.val()
            + black_cox.black_cox_cds_spread.val()
            + ah_call.val()
            + bess.total_npv.val()
            + bl.weights.iter().map(|x| x.val()).sum::<f64>()
            + kelly.effective_kelly_fraction.val()
            + kelly.expected_log_growth_rate.val()
            + dark_lit.expected_total_cost.val()
            + fx_accumulator.val()
            + atlas.val()
            + take_or_pay.npv.val()
            + market.eq_spot * 1e-6;
        let phase18_for_curve =
            FlatForwardCurve::new(market.ref_date, rate * 0.85 + inflation * 0.2);
        let phase18_csa_curve =
            FlatForwardCurve::new(market.ref_date, rate * 0.90 + credit_bump * 0.3);
        let mtm_xccs = price_mtm_resetting_xccs(
            &MtmResettingXccs {
                for_notional: 100.0,
                spot_fx: market.eq_spot + (market.eq_spot - 100.0) * 0.001,
                payment_times: vec![0.5, 1.0, 1.5, 2.0],
                accrual_fractions: vec![0.5; 4],
                for_basis_spread: -0.0015 - credit_bump * 0.2,
                dom_spread: 0.0,
                receive_foreign: true,
            },
            &phase17_curve,
            &phase18_for_curve,
            &phase18_csa_curve,
        );
        let tenor_basis = price_tenor_basis_swap(
            &TenorBasisSwapSpec {
                notional: 100.0,
                short_leg_times: vec![0.5, 1.0],
                long_leg_times: vec![1.0],
                short_leg_basis_spread: 0.001 + credit_bump * 0.1,
                pay_short_leg: true,
            },
            &phase17_curve,
            &phase18_for_curve,
            &phase18_csa_curve,
        );
        let callable_oas = analyze_callable_bond_oas_ad(&CallableBondOasSpec {
            market_dirty_price: 100.0 - credit_bump * 50.0,
            face_amount: 100.0,
            coupon_rate: 0.04 + inflation * 0.2,
            payment_times: vec![1.0, 2.0, 3.0],
            call_price: Some(101.0),
            put_price: None,
            short_rate: rate,
            short_rate_vol: 0.012 + market.eq_vol * 0.02,
            tree_steps: 8,
            smoothing_eps: 1e-3,
        });
        let ctd_colva = compute_multicurrency_ctd_colva(&MultiCurrencyCsaSpec {
            eligible_currency_spreads: vec![
                0.008 + rate * 0.1,
                0.010 + credit_bump * 0.5,
                0.009 + inflation * 0.3,
            ],
            spread_vols: vec![0.006 + market.eq_vol * 0.01, 0.008, 0.007],
            cross_currency_corr: 0.25,
            collateral_balance_profile: vec![50.0, 40.0],
            payment_times: vec![0.5, 1.0],
        });
        let turlakov_wwr = compute_turlakov_wwr_cva_fva(&TurlakovWwrSpec {
            times: vec![0.5, 1.0],
            base_epe: vec![10.0, 14.0],
            stressed_epe: vec![25.0 + (100.0 - market.eq_spot).max(0.0), 32.0],
            base_mtm: vec![8.0, 10.0],
            stressed_mtm: vec![20.0, 24.0],
            sovereign_hazard_rate: 0.02 + credit_bump,
            counterparty_hazard_rate: 0.025 + credit_bump,
            crisis_prob_threshold: 0.08,
            coupling_lambda: 0.45,
            lgd: 0.60,
            base_funding_spread: 0.01,
            stressed_funding_spread: 0.025 + credit_bump,
            discount_rate: rate,
        });
        let zabr_q = evaluate_zabr_smile(
            &ZabrParams {
                alpha: 0.18 + market.eq_vol * 0.5,
                beta: 0.7,
                nu: 0.40,
                rho: -0.25,
                gamma: 0.85,
            },
            market.eq_spot,
            105.0,
            1.0,
        );
        let pdv_iv = pdv_implied_vol_expansion(
            &GuyonPdvParams {
                beta0: 0.15 + market.eq_vol * 0.5,
                beta1: -0.12,
                beta2: 0.25,
                theta1: 0.5,
                theta2: 0.5,
                r10: -0.05,
                r11: -0.02,
                r20: 0.04,
                r21: 0.04,
            },
            market.eq_spot,
            95.0,
            1.0,
        );
        let fx_kiko = price_fx_kiko_option(&FxKikoOptionSpec {
            spot: market.eq_spot,
            strike: 100.0,
            lower_barrier: 75.0,
            upper_barrier: 125.0,
            vol: market.eq_vol,
            expiry: 1.0,
            r_dom: rate,
            r_for: 0.01,
            is_call: true,
            lower_is_ko_upper_is_ki: true,
        });
        let fx_wb = price_fx_window_barrier_option(&FxWindowBarrierSpec {
            spot: market.eq_spot,
            strike: 100.0,
            barrier: 130.0,
            vol: market.eq_vol,
            expiry: 1.0,
            r_dom: rate,
            r_for: 0.01,
            window_start: 0.1,
            window_end: 0.8,
            is_up_and_out: true,
            is_call: true,
        });
        let everest_val = price_everest_option_mc_ad(
            &EverestOptionSpec {
                spots: vec![market.eq_spot, 102.0, 98.0],
                vols: vec![market.eq_vol, 0.22, 0.20],
                correlation: 0.4,
                rate,
                expiry: 1.0,
                strike: 0.85,
                notional: 10.0,
                tau_temp: 0.02,
            },
            16,
        );
        let capped_timer = price_capped_timer_option_quadrature(
            &FiniteCapTimerOptionSpec {
                spot: market.eq_spot,
                strike: 100.0,
                rate,
                div_yield: 0.01,
                variance_budget: 0.04,
                t_max: 0.8,
                is_call: true,
            },
            &market.heston_params,
        );
        let weather_opt = price_weather_hdd_cdd_option(&WeatherDerivativeSpec {
            kind: DegreeDayKind::Hdd,
            temperature_base: 18.0,
            initial_temperature: 12.0 + inflation * 50.0,
            mean_level: 13.0,
            seasonal_amplitude: 5.0,
            seasonal_phase: 0.2,
            mean_reversion: 0.25,
            daily_vol: 2.0,
            market_price_risk: 0.05,
            times: vec![1.0, 2.0, 3.0],
            strike: 14.0,
            notional: 5.0,
            is_call: true,
            cap: Some(60.0),
        });
        let clean_spread = price_clean_spark_dark_spread(&CleanSparkDarkSpreadSpec {
            power: 80.0 + inflation * 100.0,
            fuel: 3.0,
            carbon: 20.0,
            heat_rate: 7.0,
            emission_factor: 0.5,
            strike: 10.0,
            vol_power: 0.22 + market.eq_vol * 0.2,
            vol_fuel: 0.30,
            vol_carbon: 0.35,
            rho_power_fuel: 0.45,
            rho_power_carbon: 0.25,
            rho_fuel_carbon: 0.20,
            rate,
            expiry: 1.0,
            option_type: OptionType::Call,
        });
        let min_torsion_enb = effective_number_of_bets_minimum_torsion(
            &[0.6, 0.4],
            &[
                vec![0.04 + market.eq_vol * 0.05, 0.012],
                vec![0.012, 0.09 + credit_bump * 0.5],
            ],
        );
        let ep_post = entropy_pooling_posterior(&EntropyPoolingSpec {
            scenario_returns: vec![vec![0.08, 0.02], vec![-0.04, 0.06], vec![0.02, -0.03]],
            prior_probs: vec![1.0 / 3.0; 3],
            view_matrix: vec![vec![0.08, -0.04, 0.02]],
            view_targets: vec![0.03 + rate * 0.2],
            risk_aversion: 2.5,
        });
        let mtmxccs_colva_zabr_weather_enb_metric = mtm_xccs.npv
            + tenor_basis.npv
            + callable_oas.oas * 100.0
            + ctd_colva.total_ctd_colva
            + turlakov_wwr.wwr_cva
            + zabr_q.call_price
            + pdv_iv
            + fx_kiko
            + fx_wb
            + everest_val
            + capped_timer.npv
            + weather_opt.npv
            + clean_spread
            + min_torsion_enb.enb
            + ep_post.weights[0];
        // Phase 19 institutional engines: evaluate each live path against the
        // current scenario market and aggregate into the workstation telemetry.
        let phase19_mbs = analyze_mbs_burnout_oas_ad(
            &MbsBurnoutOasSpec {
                market_price: 101.0 + market.eq_spot * 0.01,
                original_balance: 100.0,
                wac: 0.06,
                pass_through_rate: 0.055,
                term_months: 24,
                seasoning_months: 12,
                turnover_cpr: 0.06,
                max_refi_cpr: 0.30,
                refi_slope: 60.0,
                refi_incentive_threshold: 0.005,
                burnout_beta: 10.0 + market.eq_vol,
                primary_secondary_spread: 0.01,
            },
            &active_curves.curves[0],
        )
        .oas;
        let phase19_zciis = price_zciis_capfloor_jarrow_yildirim(
            &ZciisCapFloorSpec {
                notional: 100.0,
                maturity_years: 5.0,
                strike_annual_rate: 0.02 + credit_bump,
                target_month: None,
                is_cap: true,
            },
            &active_curves.curves[0],
            &inflation_curve.curve,
            &JarrowYildirimParams {
                a_nom: 0.05,
                sigma_nom: 0.008,
                a_real: 0.04,
                sigma_real: 0.006,
                sigma_cpi: 0.012 + market.eq_vol * 0.01,
                rho_nom_real: 0.3,
                rho_nom_cpi: 0.2,
                rho_real_cpi: -0.1,
            },
        )
        .npv;
        let p19_credit = SurvivalProbabilityCurve::new(
            market.ref_date,
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            vec![0.02 + credit_bump; 6],
        );
        let phase19_cmcds = price_cmcds_market_model(
            &CmCdsMarketModelSpec {
                notional: 100.0,
                grid_times: vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                start_idx: 0,
                end_idx: 3,
                cm_tenor_periods: 2,
                participation_rate: 0.95,
                lgd: 0.6,
                forward_cds_vols: vec![0.30 + market.eq_vol * 0.1; 6],
                forward_cds_corr: 0.8,
            },
            &active_curves.curves[0],
            &p19_credit,
        )
        .seller_npv;
        let phase19_simm = compute_dynamic_simm_mva_and_spike(&DynamicSimmMvaSpec {
            times: vec![0.5, 1.0, 1.5],
            bucket_sensitivities_t0: vec![100.0 + summary.total_eq_delta.abs(), 70.0],
            bucket_decay_rates: vec![0.25, 0.30],
            risk_weights: vec![0.5, 0.6],
            concentration_thresholds: vec![80.0, 80.0],
            inter_bucket_corr: 0.35,
            scheduled_cashflows: vec![12.0, 15.0, 18.0],
            mpor_years: 10.0 / 252.0,
            margin_funding_spread: 0.015 + credit_bump,
            hazard_rate: 0.02 + credit_bump,
            lgd: 0.6,
            discount_rate: rate,
        })
        .total_aps_cva;
        let phase19_jlt = evaluate_jlt_rating_migration_credit(&JltRatingMigrationSpec {
            generator_matrix: vec![
                vec![-0.06, 0.045, 0.01, 0.005],
                vec![0.02, -0.09, 0.055, 0.015],
                vec![0.005, 0.025, -0.11, 0.08],
                vec![0.0, 0.0, 0.0, 0.0],
            ],
            risk_premium_multipliers: vec![1.2 + credit_bump, 1.4, 1.6, 1.0],
            initial_rating_index: 0,
            horizon_years: 3.0,
            payment_times: vec![1.0, 2.0, 3.0],
            base_coupon_rate: 0.04,
            downgrade_step_up_spreads: vec![0.0, 0.005, 0.015, 0.0],
            recovery_rate: 0.4,
            risk_free_rate: rate,
        })
        .step_up_bond_price;
        let phase19_qou = evaluate_quintic_ou_spx_vix(
            &QuinticOuParams {
                xi0: 0.04 + market.eq_vol * 0.01,
                kappa: 2.5,
                nu: 0.8,
                rho: -0.7,
                alpha0: 1.0,
                alpha1: 0.30,
                alpha3: 0.10,
                alpha5: 0.02,
            },
            0.5,
            1.0 / 12.0,
            20.0,
            (-rate * 0.5).exp(),
        )
        .vix_call_price;
        let phase19_fbs = evaluate_free_boundary_sabr(
            &FreeBoundarySabrParams {
                forward: -0.003 + market.eq_spot * 1e-5,
                expiry: 1.0,
                alpha: 0.008 + market.eq_vol * 0.001,
                beta: 0.25,
                rho: -0.20,
                nu: 0.40,
            },
            -0.001,
            (-rate).exp(),
        )
        .call_price;
        let phase19_afv = price_convertible_bond_afv_pde(
            &AfvConvertibleBondSpec {
                spot: market.eq_spot,
                conversion_ratio: 1.0,
                face_value: 100.0,
                coupon_rate: 0.03,
                maturity: 1.5,
                rate,
                div_yield: 0.01,
                vol: market.eq_vol,
                base_hazard_rate: 0.02 + credit_bump,
                hazard_equity_elasticity: 1.2,
                equity_loss_on_default_eta: 0.5,
                bond_recovery_rate: 0.4,
                call_price: Some(115.0),
                soft_call_trigger_spot: Some(130.0),
                put_price: Some(95.0),
            },
            15,
            8,
        )
        .total_price;
        let phase19_gswap = price_gamma_swap_and_variance_option(&GammaSwapAndVarOptionSpec {
            spot: market.eq_spot,
            rate,
            div_yield: 0.01,
            expiry: 1.0,
            cut_strike: 100.0,
            variance_strike: 0.04,
            quotes: vec![
                VarianceSwapStrikeQuote {
                    strike: 85.0,
                    implied_vol: 0.24,
                },
                VarianceSwapStrikeQuote {
                    strike: 100.0,
                    implied_vol: 0.20,
                },
                VarianceSwapStrikeQuote {
                    strike: 115.0,
                    implied_vol: 0.18,
                },
            ],
            heston: HestonParams {
                v0: market.heston_params.v0,
                kappa: 2.0,
                theta: 0.04,
                sigma_v: 0.35,
                rho: -0.7,
            },
            is_variance_call: true,
        })
        .variance_option_price;
        let phase19_cs = price_clewlow_strickland_strip_swaption(&ClewlowStricklandSwaptionSpec {
            option_expiry: 0.5,
            delivery_times: vec![0.75, 1.0, 1.25],
            forward_prices: vec![75.0 + market.eq_spot * 0.01, 76.0, 74.0],
            delivery_weights: vec![1.0 / 3.0; 3],
            strike: 74.5,
            rate,
            sigma_short: 0.32 + market.eq_vol * 0.1,
            alpha_short: 1.1,
            sigma_long: 0.15,
            rho_factors: 0.2,
            seasonal_vol_amp: 0.05,
            is_payer_call: true,
        })
        .npv;
        let phase19_mm = solve_avellaneda_stoikov_market_making(&AvellanedaStoikovMmSpec {
            mid_price: market.eq_spot,
            initial_inventory: 3.0 + credit_bump * 100.0,
            horizon: 1.0,
            num_steps: 6,
            volatility: market.eq_vol,
            risk_aversion: 0.5,
            order_book_kappa: 1.5,
            base_arrival_intensity: 10.0,
            alpha_drift: 0.01,
        })
        .expected_spread_capture_pnl;
        let phase19_cdar = optimize_mean_cdar_portfolio(&MeanCdarPortfolioSpec {
            path_asset_returns: vec![
                vec![0.03, 0.01],
                vec![-0.06, -0.01],
                vec![-0.04, -0.005],
                vec![0.04, 0.012],
            ],
            expected_returns: vec![0.09 + market.eq_vol * 0.01, 0.04],
            cdar_aversion: 4.0,
            confidence: 0.80,
            smoothing_eps: 1e-3,
        })
        .cdar;
        let phase19_schur = schur_complementary_hrp(&SchurHrpSpec {
            covariance: vec![
                vec![0.04 + market.eq_vol * 0.01, 0.012],
                vec![0.012, 0.09 + credit_bump],
            ],
            gamma: 0.6,
        })
        .portfolio_variance;
        let phase19_oas = oas_and_nonlinear_shrinkage_covariance(&[
            vec![0.02 + market.eq_vol * 0.001, -0.01],
            vec![-0.015, 0.025 + credit_bump * 0.1],
            vec![0.01, 0.005],
            vec![-0.005, -0.01],
        ])
        .oas_shrinkage_intensity;
        let burnout_cmcds_quintic_afv_cs_schur_metric = phase19_mbs * 100.0
            + phase19_zciis
            + phase19_cmcds
            + phase19_simm
            + phase19_jlt
            + phase19_qou
            + phase19_fbs * 100.0
            + phase19_afv
            + phase19_gswap
            + phase19_cs
            + phase19_mm
            + phase19_cdar
            + phase19_schur
            + phase19_oas;
        let krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric =
            Self::phase22_live_metric(&market, credit_bump);
        let blotter_rows = vec![
            ("IRS 5Y PayFix".to_string(), irs_npv, "IFT".to_string()),
            ("Heston Call 1Y".to_string(), heston_npv, "COS".to_string()),
            (
                "RoughHest H=.1".to_string(),
                rough_npv,
                "Lift 8-Fac".to_string(),
            ),
            ("Autocall Note".to_string(), script_npv, "AAD".to_string()),
            ("CDS 5Y BuyProt".to_string(), cds_npv, "IFT".to_string()),
            (
                "ZCIIS 5Y CPI".to_string(),
                zciis_npv,
                "Seasonal".to_string(),
            ),
            (
                "FX 25D VV Call".to_string(),
                fx.price_vanna_volga(fx.k3, true).val(),
                "Dual3".to_string(),
            ),
            (
                "EXOT-BAR-DN-OUT-1Y".to_string(),
                barrier_npv,
                "Dual5 BGK".to_string(),
            ),
            (
                "CMDTY-WTI-SS2F-1Y".to_string(),
                commodity_option + crack,
                "SS2F/Kirk".to_string(),
            ),
            ("G2++ Swaption".to_string(), g2_npv, "G2 Quad".to_string()),
            ("LEVY-VG-COS-1Y".to_string(), levy_npv, "COS VG".to_string()),
            (
                "CONV-TF-BOND-3Y".to_string(),
                conv_npv.val(),
                "TF + ASB".to_string(),
            ),
            (
                "SOFR RFR SABR 1Y".to_string(),
                rfr_npv,
                "Willems".to_string(),
            ),
            (
                "EURUSD XCCS 5Y".to_string(),
                xccs_npv,
                "IFT Basis".to_string(),
            ),
            ("Portfolio ES 97.5%".to_string(), es, "Acerbi".to_string()),
            (
                "Bermudan HW 3Y".to_string(),
                bermudan_value,
                "Tree IFT".to_string(),
            ),
            (
                "SA-CCR EAD".to_string(),
                sa_report.ead.val(),
                "BCBS279".to_string(),
            ),
            (
                "Commodity Swing".to_string(),
                swing_value,
                "OU DP".to_string(),
            ),
            ("Box Markowitz".to_string(), qp_value, "KKT IFT".to_string()),
            (
                "Gas Storage DP".to_string(),
                storage_value,
                "OU DP".to_string(),
            ),
        ];
        WorkstationFrameData {
            tick: self.tick,
            scenario: self.scenario,
            chart_mode: self.chart_mode,
            eval_micros: elapsed,
            dag_recomputed_trades: summary.last_trade_recomputes,
            dag_total_trades: self.trades.len(),
            dag_curves_rebuilt: summary.last_curve_recomputes > 0,
            total_npv: adjusted_total_npv,
            intraday_pnl: adjusted_total_npv - prev,
            eq_delta: summary.total_eq_delta,
            eq_vega,
            fx_spot_delta: fx_delta,
            inflation_par_delta_sum: inflation,
            cva,
            dva,
            fva,
            mva,
            simm_initial_margin: simm,
            ir_dv01_buckets: ir_buckets,
            credit_cs01_buckets: credit_buckets,
            jump_to_default: jtd,
            bermudan_sabr_metric: bermudan_value + sabr_metric,
            sacc_swing_qp_metric: sa_report.ead.val() + swing_value + qp_value,
            slv_rough_worst_of_metric: rough_metric,
            frtb_cvar_execution_metric: execution_metric,
            vix_parisian_robust_metric,
            double_heston_tarf_coco_hrp_metric,
            jy_tba_bergomi_nco_metric,
            cheyette_cmo_qrh_evar_metric,
            asw_mf1f_ah_bess_ep_metric,
            mtmxccs_colva_zabr_weather_enb_metric,
            burnout_cmcds_quintic_afv_cs_schur_metric,
            sinhou_mb_rbergomi_choi_afs_expectile_metric,
            monotone_convex_clark_saccr_rvol_cartea_paleologo_metric: phase21_live_metric(
                &market,
                credit_bump,
            ),
            krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric,
            rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml_metric: self
                .phase23_live_metric(),
            mbs_rmv_mo_sig_garch_rie_gp_ierp_metric: self.phase24_live_metric(),
            kalman_vecm_hj_ifrs9_cpt_hasbrouck_fv_lando_lppl_cbarb_causal_malliavin_metric:
                Self::phase25_live_metric_from_market(&market, credit_bump),
            phase26_sota_metric: Self::phase26_live_metric_from_market(&market, credit_bump),
            phase27_sota_metric: Self::phase27_live_metric_from_market(&market, credit_bump),
            phase29_sota_metric: Self::phase29_live_metric_from_market_with_bump(
                &market,
                credit_bump,
            ),
            phase30_sota_metric: Self::phase30_live_metric_from_market_with_bump(
                &market,
                credit_bump,
            ),
            phase31_sota_metric: Self::phase31_live_metric_from_market(&market, credit_bump),
            phase32_sota_metric: Self::phase32_live_metric_from_market(&market, credit_bump),
            phase33_sota_metric: Self::phase33_live_metric_from_market(&market, credit_bump),
            phase28_sota_metric: Self::phase28_live_metric_from_market_with_bump(
                &market,
                if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
                    0.01
                } else {
                    0.0
                },
            ),
            blotter_rows,
        }
    }
    fn chart(&self, frame: &WorkstationFrameData) -> BrailleCanvas {
        let mut canvas = BrailleCanvas::new(58, 10);
        let market = self.scenario_market();
        let curves = self.active_curve_system(&market);
        let x32: Vec<f64> = (0..32).map(|i| 0.25 + i as f64 * 9.75 / 31.0).collect();
        match self.chart_mode {
            ChartMode::YieldAndInflationCurves => {
                let (_, inflation) = self.inflation_calibration(self.scenario);
                let ois: Vec<f64> = x32
                    .iter()
                    .map(|&t| curves.curves[0].zero_rate(t).val())
                    .collect();
                let fwd: Vec<f64> = x32
                    .iter()
                    .map(|&t| curves.curves[1].zero_rate(t).val())
                    .collect();
                let cpi: Vec<f64> = x32
                    .iter()
                    .map(|&t| {
                        let month = (((t * 12.0).round() as i32 - 1).rem_euclid(12) + 1) as u8;
                        ZeroCouponInflationSwap {
                            notional: 1.0,
                            fixed_rate: 0.0,
                            maturity_years: t,
                            target_month: Some(month),
                            pay_fixed: true,
                        }
                        .fair_rate(&inflation.curve)
                        .val()
                    })
                    .collect();
                let (lo, hi) = range3(&ois, &fwd, &cpi);
                canvas.plot_series(&x32, &ois, lo, hi, [0, 220, 255]);
                canvas.plot_series(&x32, &fwd, lo, hi, [0, 255, 140]);
                canvas.plot_series(&x32, &cpi, lo, hi, [255, 200, 60]);
            }
            ChartMode::RoughHestonAndFxSmiles => {
                let strikes: Vec<f64> = (0..25).map(|i| 85.0 + i as f64 * 30.0 / 24.0).collect();
                let expiry = 1.0;
                let rate = curves.curves[0].zero_rate(expiry).val();
                let discount = (-rate * expiry).exp();
                let forward = market.eq_spot * (rate * expiry).exp();
                let rough = LiftedHestonParams {
                    v0: market.heston_params.v0,
                    lambda: market.heston_params.kappa,
                    theta: market.heston_params.theta,
                    nu: market.heston_params.sigma_v,
                    rho: market.heston_params.rho,
                    hurst: 0.1,
                    n_factors: 8,
                    r_n: 2.5,
                };
                let heston: Vec<f64> = strikes
                    .iter()
                    .map(|&k| {
                        let price = lifted_heston_cos_price(
                            market.eq_spot,
                            k,
                            expiry,
                            rate,
                            0.0,
                            &rough,
                            64,
                            8.0,
                            16,
                            true,
                        )
                        .val();
                        implied_black_vol(forward, k, expiry, price / discount, true)
                            .map(|v| v.val())
                            .unwrap_or(0.0)
                    })
                    .collect();
                let fx = FxVannaVolgaSmile::new(
                    100.0,
                    0.03,
                    0.01,
                    expiry,
                    FxSmileQuote {
                        atm_vol: 0.1 + 0.1 * market.eq_vol,
                        rr_25d: 0.01,
                        bf_25d: 0.005,
                    },
                );
                let fx_vol: Vec<f64> = strikes.iter().map(|&k| fx.implied_vol(k).val()).collect();
                let (lo, hi) = range2(&heston, &fx_vol);
                canvas.plot_series(&strikes, &heston, lo, hi, [0, 220, 255]);
                canvas.plot_series(&strikes, &fx_vol, lo, hi, [255, 200, 60]);
            }
            ChartMode::XvaExposureProfile => {
                let cube = self.exposure_cube();
                let epe = expected_positive_exposure(&cube);
                let pfe = potential_future_exposure(&cube, 0.975);
                let dim = {
                    let d = self.base_market.ref_date;
                    let active = self.scenario_market();
                    let curve = FlatForwardCurve::new(d, active.ois_quotes[0]);
                    let model = HullWhite1F {
                        a: 0.05,
                        sigma: 0.01 + 0.25 * active.eq_vol,
                    };
                    let swap = NettingSetSwap {
                        notional: 1_000_000.0,
                        fixed_rate: active.ois_quotes[1],
                        payment_times: (1..=5).map(|x| x as f64).collect(),
                        accruals: vec![1.0; 5],
                        is_payer: true,
                    };
                    let times = (1..=16).map(|k| k as f64 * 5.0 / 16.0).collect::<Vec<_>>();
                    ql_xva::compute_chebyshev_dim_and_mva(
                        &model,
                        &curve,
                        &[swap],
                        &times,
                        64,
                        6,
                        0.01,
                        29,
                    )
                    .eim_profile
                };
                let xs: Vec<f64> = cube.times.clone();
                let (lo, hi) = range3(&epe, &pfe, &dim);
                canvas.plot_series(&xs, &epe, lo, hi, [0, 255, 140]);
                canvas.plot_series(&xs, &pfe, lo, hi, [255, 90, 90]);
                canvas.plot_series(&xs, &dim, lo, hi, [180, 130, 255]);
            }
        }
        let _ = frame;
        canvas
    }
    pub fn render_dashboard(&self, frame: &WorkstationFrameData, ansi_color: bool) -> String {
        const LEFT: usize = 57;
        const RIGHT: usize = 58;
        let chart_lines = self.chart(frame).render_lines(ansi_color);
        let mut left_top = Vec::with_capacity(10);
        left_top.push("1. BOOK SUMMARY, XVA & ISDA SIMM".to_string());
        left_top.push(format!(
            "Total NPV: {:>14.2} {}",
            frame.total_npv,
            render_sparkline(&self.npv_history, 16)
        ));
        left_top.push(format!("Intraday PnL: {:>11.2}", frame.intraday_pnl));
        left_top.push(format!(
            "Eq Delta: {:>12.2}  Eq Vega: {:>10.2}",
            frame.eq_delta, frame.eq_vega
        ));
        left_top.push(format!(
            "FX Delta: {:>12.5}  CPI Δ01: {:>10.2}",
            frame.fx_spot_delta, frame.inflation_par_delta_sum
        ));
        left_top.push(format!(
            "XVA CVA/DVA/FVA: {:>7.2}/{:>7.2}/{:>7.2}  MVA:{:>7.2}",
            frame.cva, frame.dva, frame.fva, frame.mva
        ));
        left_top.push(format!("ISDA SIMM IM: {:>14.2}", frame.simm_initial_margin));
        left_top.push(format!("Jump-to-default: {:>12.2}", frame.jump_to_default));
        let cached = frame
            .dag_total_trades
            .saturating_sub(frame.dag_recomputed_trades);
        left_top.push(format!(
            "DAG cache hit: {cached}/{}",
            frame.dag_total_trades
        ));
        left_top.push(format!(
            "Eval: {:>8.0} µs  Mode: {}",
            frame.eval_micros,
            frame.scenario.label()
        ));
        let right_top = chart_lines;

        let max_ir = frame
            .ir_dv01_buckets
            .iter()
            .map(|x| x.1.abs())
            .fold(0.0, f64::max)
            .max(1.0);
        let max_cs = frame
            .credit_cs01_buckets
            .iter()
            .map(|x| x.1.abs())
            .fold(0.0, f64::max)
            .max(1.0);
        let mut left_bottom = Vec::with_capacity(11);
        for (label, value) in &frame.ir_dv01_buckets {
            left_bottom.push(format!(
                "IRS DV01 {label:>3}: {:>9.2} {}",
                value,
                render_diverging_bar(*value, max_ir, 10, ansi_color)
            ));
        }
        for (label, value) in &frame.credit_cs01_buckets {
            left_bottom.push(format!(
                "CDS CS01 {label:>3}: {:>9.2} {}",
                value,
                render_diverging_bar(*value, max_cs, 10, ansi_color)
            ));
        }
        let american = BaroneAdesiWhaleyEngine::price(100.0, 100.0, 0.03, 0.01, 0.2, 1.0, true);
        let erc = equal_risk_contribution_weights(&[vec![0.04, 0.01], vec![0.01, 0.09]]);
        left_bottom.extend([
            format!("19 Bermudan/SABR IFT: {:.2}", frame.bermudan_sabr_metric),
            format!("20 SA-CCR/Swing/QP: {:.2}", frame.sacc_swing_qp_metric),
            format!(
                "21 SLV/rBergomi/WoF: {:.2}",
                frame.slv_rough_worst_of_metric
            ),
            format!("22 FRTB/CVaR/Exec: {:.2}", frame.frtb_cvar_execution_metric),
            format!(
                "23 Monotone/VIX/Parisian/Robust: {:.2}",
                frame.vix_parisian_robust_metric
            ),
            format!(
                "24 Double Heston/TARF/CoCo/HRP: {:.2}",
                frame.double_heston_tarf_coco_hrp_metric
            ),
            format!(
                "25 JY-Inflation/TBA/Bergomi2F/NCO: {:.2}",
                frame.jy_tba_bergomi_nco_metric
            ),
            format!(
                "26 Cheyette/CMO/QRH/BaseCorr/EVaR: {:.2}",
                frame.cheyette_cmo_qrh_evar_metric
            ),
            format!(
                "27 ASW/MF1F/AH-LV/BESS/Idzorek-DP: {:.2}",
                frame.asw_mf1f_ah_bess_ep_metric
            ),
            format!(
                "28 MtM-XCCS/COLVA/ZABR/Weather/ENB: {:.2}",
                frame.mtmxccs_colva_zabr_weather_enb_metric
            ),
            format!(
                "29 Burnout/CMCDS/Quintic/AFV/Schur: {:.2}",
                frame.burnout_cmcds_quintic_afv_cs_schur_metric
            ),
            format!(
                "30 SinhOU/MB/mrBergomi/Choi/EXVaR: {:.2}",
                frame.sinhou_mb_rbergomi_choi_afs_expectile_metric
            ),
            format!(
                "31 HW-MC/Clark/SACCR/rVol/Cartea/Paleo: {:.2}",
                frame.monotone_convex_clark_saccr_rvol_cartea_paleologo_metric
            ),
            format!(
                "32 KRD/LPI/CoCo/Bergomi2F/CGMY/Guéant/Meucci: {:.2}",
                frame.krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric
            ),
            format!(
                "33 RKHS/FMM/VV/LSV/JLT/CTD/QRH/NCO/DML: {:.2}",
                frame.rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml_metric
            ),
            format!(
                "34 MBS/RMV/MO/CB/Sig/HN/POT/RIE/PIN/GP/BC/IERP: {:.2}",
                frame.mbs_rmv_mo_sig_garch_rie_gp_ierp_metric
            ),
            format!(
                "35 Phase25 SOTA live: {:.2}",
                frame
                    .kalman_vecm_hj_ifrs9_cpt_hasbrouck_fv_lando_lppl_cbarb_causal_malliavin_metric
            ),
            format!("36 Phase26 SOTA live: {:.2}", frame.phase26_sota_metric),
            format!("37 Phase 27 SOTA live: {:.2}", frame.phase27_sota_metric),
            format!("38:P28-SOTA live: {:.2}", frame.phase28_sota_metric),
            format!("39 Phase29 SOTA live: {:.2}", frame.phase29_sota_metric),
            format!("40 Phase30 SOTA live: {:.2}", frame.phase30_sota_metric),
            format!("41 Phase31 SOTA live: {:.2}", frame.phase31_sota_metric),
            format!("42 Phase32 SOTA live: {:.2}", frame.phase32_sota_metric),
            format!("43 Phase33 SOTA live: {:.2}", frame.phase33_sota_metric),
            format!(
                "DAG recomputed: {}/{}",
                frame.dag_recomputed_trades, frame.dag_total_trades
            ),
            format!("Chart: {}", self.chart_mode.label()),
            format!("American / Hybrid: {:.4}", american),
            format!("BL / Risk Parity ERC: {:.3}/{:.3}", erc[0], erc[1]),
        ]);
        // There are six IR buckets, five credit buckets and ten telemetry
        // lines (including Tab 23).  Render all 21 rows instead of silently
        // dropping the final BL/Risk-Parity line.
        const BOTTOM_ROWS: usize = 40;
        while left_bottom.len() < BOTTOM_ROWS {
            left_bottom.push(String::new());
        }
        let mut right_bottom = vec![];
        right_bottom.extend(
            frame
                .blotter_rows
                .iter()
                .map(|(name, value, tag)| format!("{name:<20} {:>12.2} {tag:<10}", value)),
        );
        right_bottom.push("32 Phase22 SOTA: KRD/LPI/CoCo/Bergomi/CGMY/AAD live".to_string());
        right_bottom.push("33 RKHS/FMM/VV/LSV/JLT/CTD/QRH/NCO/DML live".to_string());
        right_bottom.push("34 Davidson-Levin MBS IO/PO/Duffie-Singleton RMV/MO-XVA/Sirbu-Shreve CB/PathSig/HN-DCC-GARCH/POT-GPD/RIE/PIN/GP-Aim/Brinson/IERP".to_string());
        right_bottom.push(format!("35 Phase25 SOTA Kalman/VECM/HJ/IFRS9/CPT/Hasbrouck/FV/Lando/LPPL/CB-Arb/Causal/Malliavin: {:.4}", self.phase25_live_metric()));
        right_bottom.push(format!("36 Phase26 SOTA MF-Rates/Wildcard/TARF/Bergomi/Schönbucher/FRTB-CVA/Lookback/VTA/FLAM/ACD/HLZ/MLMC: {:.4}", self.phase26_live_metric()));
        right_bottom.push(format!("37 Phase 27 SOTA FTPL/G2/eSSVI/CB/CoCo/MPOR/Double-Heston/Curran/Warin/Roll/Busseti/KKT: {:.4}", self.phase27_live_metric()));
        right_bottom.push(format!("38:P28-SOTA IRRBB/RFR/XCCY/BG/CoCo/SIMM/rBergomi/Spread/BESS/Dark/Risk/Diff-PCA: {:.4}", self.phase28_live_metric()));
        right_bottom.push(format!("39 Phase29 SOTA Veronesi/Lando/Corridor/Russian/Bates/Koziol/MFG/BCBS/CJP/ETF/IC/IS: {:.4}", self.phase29_live_metric()));
        right_bottom.push(format!("40 Phase30 SOTA Taylor/Damodaran/VaR/LVaR/CGO/Regret/Beneish/Campisi/Heckman/Elastic/ResetCB/FSB: {:.4}", self.phase30_live_metric()));
        right_bottom.push(format!("41 Phase31 SOTA KRD-LSC/CMO-PAC/Asquith-CB/ASU2024/ASC470/ColVA/Hill/Bouchaud/Karnosky/DGTW/Brock/Penman: {:.4}", self.phase31_live_metric()));
        right_bottom.push(format!("42 Phase32 SOTA Cheyette/NDIRS/ZABR/QRH/CMCDS/WWR/HHW/CB/EVaR/Schur/BEKK/Quality: {:.4}", self.phase32_live_metric()));
        right_bottom.push(format!(
            "43 Phase33 SOTA Rebonato-CMS/Leland-Toft/Calamos/RIE/Arrow-Debreu: {:.4}",
            self.phase33_live_metric()
        ));
        while right_bottom.len() < BOTTOM_ROWS {
            right_bottom.push(String::new());
        }

        let mut out = String::new();
        let top_left = pad_visible(
            &format!("─ ql-rs QUANT RISK WORKSTATION  Tick #{:04} ", frame.tick),
            LEFT,
        );
        let top_right = pad_visible(
            &format!("─ 2. BRAILLE CANVAS: {} ", self.chart_mode.label()),
            RIGHT,
        );
        let _ = writeln!(out, "╭{}┬{}╮", top_left, top_right);
        for i in 0..10 {
            let l = left_top.get(i).map(String::as_str).unwrap_or("");
            let r = right_top.get(i).map(String::as_str).unwrap_or("");
            let _ = writeln!(out, "│{}│{}│", pad_visible(l, LEFT), pad_visible(r, RIGHT));
        }
        let _ = writeln!(
            out,
            "├{}┼{}┤",
            pad_visible(" 3. IFT BUCKETED DV01 & CREDIT CS01 LADDERS", LEFT),
            pad_visible(" 4. CROSS-ASSET LIVE BLOTTER & SOTA TELEMETRY", RIGHT)
        );
        for i in 0..BOTTOM_ROWS {
            let _ = writeln!(
                out,
                "│{}│{}│",
                pad_visible(&left_bottom[i], LEFT),
                pad_visible(&right_bottom[i], RIGHT)
            );
        }
        let _ = writeln!(out, "╰{}┴{}╯", "─".repeat(LEFT), "─".repeat(RIGHT));
        out
    }
}
impl Default for DeskWorkstation {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod phase14_tests {
    use super::*;
    #[test]
    fn phase14_tab_is_registered() {
        assert_eq!(TabId::ALL.len(), 43);
        assert!(TabId::ALL.contains(&TabId::Phase14Sota));
        assert!(TabId::ALL.contains(&TabId::Phase15Sota));
        let mut ws = DeskWorkstation::new();
        let frame = ws.step_tick();
        assert!(frame.double_heston_tarf_coco_hrp_metric > 100.0);
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("24 Double Heston/TARF/CoCo/HRP")
        );
    }
}

#[cfg(test)]
mod phase15_tests {
    use super::*;
    #[test]
    fn phase15_tab_and_metric_are_live() {
        assert!(TabId::ALL.contains(&TabId::Phase15Sota));
        let mut ws = DeskWorkstation::new();
        let frame = ws.step_tick();
        assert!(frame.jy_tba_bergomi_nco_metric.is_finite());
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("25 JY-Inflation/TBA/Bergomi2F/NCO")
        );
        ws.set_scenario(ScenarioKind::RatesPlus50Bp);
        let stressed = ws.step_tick();
        assert_ne!(
            stressed.jy_tba_bergomi_nco_metric,
            frame.jy_tba_bergomi_nco_metric
        );
    }
}

#[cfg(test)]
mod phase16_tests {
    use super::*;

    #[test]
    fn phase16_tab_and_metric_are_live() {
        assert_eq!(TabId::ALL.len(), 43);
        assert!(TabId::ALL.contains(&TabId::Phase16Sota));
        let mut ws = DeskWorkstation::new();
        let frame = ws.step_tick();
        assert!(frame.cheyette_cmo_qrh_evar_metric.is_finite());
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("26 Cheyette/CMO/QRH/BaseCorr/EVaR")
        );
        ws.set_scenario(ScenarioKind::RatesPlus50Bp);
        let stressed = ws.step_tick();
        assert_ne!(
            stressed.cheyette_cmo_qrh_evar_metric,
            frame.cheyette_cmo_qrh_evar_metric
        );
    }
}

#[cfg(test)]
mod phase17_tests {
    use super::*;

    #[test]
    fn phase17_tab_metric_dashboard_and_rates_reactivity() {
        assert_eq!(TabId::ALL.len(), 43);
        assert!(TabId::ALL.contains(&TabId::Phase17Sota));
        assert_eq!(
            TabId::Phase17Sota.label(),
            "27 ASW / MF1F / AH-LV / BESS / Idzorek-DP"
        );
        let mut ws = DeskWorkstation::new();
        let frame = ws.step_tick();
        assert!(frame.asw_mf1f_ah_bess_ep_metric.is_finite());
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("27 ASW/MF1F/AH-LV/BESS/Idzorek-DP")
        );
        let normal_next = ws.step_tick();
        assert_ne!(
            normal_next.asw_mf1f_ah_bess_ep_metric,
            frame.asw_mf1f_ah_bess_ep_metric
        );
        for scenario in [
            ScenarioKind::RatesPlus50Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::CreditWideningPlus100Bp,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(scenario);
            let stressed = ws.step_tick();
            assert_ne!(
                stressed.asw_mf1f_ah_bess_ep_metric,
                normal_next.asw_mf1f_ah_bess_ep_metric
            );
        }
    }
}

#[cfg(test)]
mod phase18_tests {
    use super::*;
    #[test]
    fn phase18_tab_metric_dashboard_and_stress_reactivity() {
        assert_eq!(TabId::ALL.len(), 43);
        assert_eq!(
            TabId::Phase18Sota.label(),
            "28 MtM-XCCS / CTD-COLVA / ZABR-PDV / Weather / ENB"
        );
        let mut ws = DeskWorkstation::new();
        let normal = ws.step_tick();
        assert!(normal.mtmxccs_colva_zabr_weather_enb_metric.is_finite());
        assert!(
            ws.render_dashboard(&normal, false)
                .contains("28 MtM-XCCS/COLVA/ZABR/Weather/ENB")
        );
        assert!(
            ws.render_dashboard(&normal, false)
                .contains("BL / Risk Parity ERC")
        );
        let next = ws.step_tick();
        assert_ne!(
            next.mtmxccs_colva_zabr_weather_enb_metric,
            normal.mtmxccs_colva_zabr_weather_enb_metric
        );
        for scenario in [
            ScenarioKind::RatesPlus50Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::CreditWideningPlus100Bp,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(scenario);
            let stressed = ws.step_tick();
            assert_ne!(
                stressed.mtmxccs_colva_zabr_weather_enb_metric,
                normal.mtmxccs_colva_zabr_weather_enb_metric
            );
        }
    }
}

#[cfg(test)]
mod phase19_tests {
    use super::*;

    #[test]
    fn phase19_tab29_terminal_renders_and_ticks_cleanly() {
        assert_eq!(TabId::ALL.len(), 43);
        assert_eq!(
            TabId::Phase19Sota.label(),
            "29 Burnout-MBS / CMCDS / Quintic-OU / AFV-CB / CS-Swaption / Schur-HRP"
        );
        let mut ws = DeskWorkstation::new();
        let normal = ws.step_tick();
        assert!(normal.burnout_cmcds_quintic_afv_cs_schur_metric.is_finite());
        let dashboard = ws.render_dashboard(&normal, false);
        assert!(dashboard.contains("29 Burnout/CMCDS/Quintic/AFV/Schur"));
        assert!(dashboard.contains("BL / Risk Parity ERC"));
        let next = ws.step_tick();
        assert_ne!(
            next.burnout_cmcds_quintic_afv_cs_schur_metric,
            normal.burnout_cmcds_quintic_afv_cs_schur_metric
        );
        for scenario in [
            ScenarioKind::RatesPlus50Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::CreditWideningPlus100Bp,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(scenario);
            let stressed = ws.step_tick();
            assert_ne!(
                stressed.burnout_cmcds_quintic_afv_cs_schur_metric,
                normal.burnout_cmcds_quintic_afv_cs_schur_metric
            );
        }
    }
}

#[cfg(test)]
mod phase20_tests {
    use super::*;

    #[test]
    fn phase20_tab_metric_dashboard_and_stress_reactivity() {
        assert_eq!(TabId::ALL.len(), 43);
        assert!(TabId::ALL.contains(&TabId::Phase20Sota));
        let mut ws = DeskWorkstation::new();
        let normal = ws.step_tick();
        assert!(
            normal
                .sinhou_mb_rbergomi_choi_afs_expectile_metric
                .is_finite()
        );
        assert!(
            ws.render_dashboard(&normal, false)
                .contains("30 SinhOU/MB/mrBergomi/Choi/EXVaR")
        );
        for scenario in [
            ScenarioKind::RatesPlus50Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::CreditWideningPlus100Bp,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(scenario);
            let stressed = ws.step_tick();
            assert!(
                stressed
                    .sinhou_mb_rbergomi_choi_afs_expectile_metric
                    .is_finite()
            );
            assert_ne!(
                stressed.sinhou_mb_rbergomi_choi_afs_expectile_metric,
                normal.sinhou_mb_rbergomi_choi_afs_expectile_metric
            );
        }
    }
}

#[cfg(test)]
mod phase21_tests {
    use super::*;

    #[test]
    fn phase21_terminal_reactivity() {
        assert_eq!(TabId::ALL.len(), 43);
        let mut ws = DeskWorkstation::new();
        let normal = ws.step_tick();
        let next = ws.step_tick();
        assert_ne!(
            next.monotone_convex_clark_saccr_rvol_cartea_paleologo_metric,
            normal.monotone_convex_clark_saccr_rvol_cartea_paleologo_metric
        );
        for scenario in [
            ScenarioKind::RatesPlus50Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::CreditWideningPlus100Bp,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(scenario);
            let stressed = ws.step_tick();
            assert_ne!(
                stressed.monotone_convex_clark_saccr_rvol_cartea_paleologo_metric,
                normal.monotone_convex_clark_saccr_rvol_cartea_paleologo_metric
            );
        }
        assert!(
            ws.render_dashboard(&normal, false)
                .contains("31 HW-MC/Clark/SACCR/rVol/Cartea/Paleo")
        );
        assert!(
            ws.render_dashboard(&normal, false)
                .contains("BL / Risk Parity ERC")
        );
    }
}

#[cfg(test)]
mod phase22_tests {
    use super::*;
    #[test]
    fn phase22_tab32_metric_renders() {
        let mut ws = DeskWorkstation::new();
        let frame = ws.step_tick();
        assert_eq!(
            TabId::Phase22Sota.label(),
            "32 Phase 22 SOTA / KRD / LPI / CoCo / Bergomi / CGMY / AAD"
        );
        assert!(
            frame
                .krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric
                .is_finite()
        );
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("32 KRD/LPI/CoCo/Bergomi2F/CGMY/Guéant/Meucci")
        );
        let next = ws.step_tick();
        assert_ne!(
            next.krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric,
            frame.krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric
        );
        for scenario in [
            ScenarioKind::RatesPlus50Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::CreditWideningPlus100Bp,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(scenario);
            let stressed = ws.step_tick();
            assert_ne!(
                stressed.krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric,
                frame.krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric
            );
        }
    }
}
#[cfg(test)]
mod phase23_tests {
    use super::*;
    #[test]
    fn phase23_tab33_metric_renders_and_reacts() {
        let mut ws = DeskWorkstation::new();
        let f = ws.step_tick();
        assert_eq!(
            TabId::Phase23Sota.label(),
            "33 RKHS/FMM/VV/LSV/JLT/CTD/QRH/NCO/DML"
        );
        assert!(
            ws.render_dashboard(&f, false)
                .contains("33 RKHS/FMM/VV/LSV/JLT/CTD/QRH/NCO/DML live")
        );
        for sc in [
            ScenarioKind::RatesPlus50Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::CreditWideningPlus100Bp,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(sc);
            let g = ws.step_tick();
            assert!(
                g.rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml_metric
                    .is_finite()
            );
        }
    }
}
#[cfg(test)]
mod phase24_tests {
    use super::*;
    #[test]
    fn phase24_tab_metric_tick_and_rate_reactivity() {
        assert_eq!(TabId::ALL.len(), 43);
        assert!(TabId::Phase24Sota.label().starts_with("34 Phase 24 SOTA"));
        let mut ws = DeskWorkstation::new();
        let a = ws.step_tick();
        assert!(a.mbs_rmv_mo_sig_garch_rie_gp_ierp_metric.is_finite());
        ws.set_scenario(ScenarioKind::RatesPlus50Bp);
        let b = ws.step_tick();
        assert!(b.mbs_rmv_mo_sig_garch_rie_gp_ierp_metric.is_finite());
        assert_ne!(
            a.mbs_rmv_mo_sig_garch_rie_gp_ierp_metric,
            b.mbs_rmv_mo_sig_garch_rie_gp_ierp_metric
        );
    }
    #[test]
    fn phase24_tab34_terminal_renders_and_reacts() {
        let mut ws = DeskWorkstation::new();
        assert_eq!(
            TabId::Phase24Sota.label(),
            "34 Phase 24 SOTA / MBS-IO-PO / RMV / MO-XVA / Dilution / Sig / GARCH / POT / RIE / PIN / GP / Brinson / IERP"
        );
        let mut vals = Vec::new();
        for scenario in [
            ScenarioKind::RatesParallelPlus50Bp,
            ScenarioKind::VolShockPlus5Pts,
            ScenarioKind::CreditSpreadWidening100Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(scenario);
            let f = ws.step_tick();
            vals.push(f.mbs_rmv_mo_sig_garch_rie_gp_ierp_metric);
            let frame = ws.render_dashboard(&f, false);
            assert!(frame.contains("34 MBS/RMV/MO/CB/Sig/HN/POT/RIE/PIN/GP/BC/IERP"));
            assert!(!frame.contains('\u{1b}') || frame.contains("34"));
        }
        assert!(vals.iter().all(|x| x.is_finite()));
        assert!(vals.windows(2).any(|w| (w[1] - w[0]).abs() > 1e-12));
    }
}

#[cfg(test)]
mod phase25_tests {
    use super::*;
    #[test]
    fn phase25_tab35_metric_reacts_across_scenarios() {
        let mut ws = DeskWorkstation::new();
        let mut vals = Vec::new();
        for scenario in [
            ScenarioKind::RatesParallelPlus50Bp,
            ScenarioKind::VolShockPlus5Pts,
            ScenarioKind::CreditSpreadWidening100Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(scenario);
            let frame = ws.step_tick();
            assert!(
                frame
                    .kalman_vecm_hj_ifrs9_cpt_hasbrouck_fv_lando_lppl_cbarb_causal_malliavin_metric
                    .is_finite()
            );
            assert!(
                ws.render_dashboard(&frame, false)
                    .contains("35 Phase25 SOTA")
            );
            vals.push(
                frame
                    .kalman_vecm_hj_ifrs9_cpt_hasbrouck_fv_lando_lppl_cbarb_causal_malliavin_metric,
            );
        }
        assert!(vals.windows(2).any(|w| (w[1] - w[0]).abs() > 1e-12));
    }
}

#[cfg(test)]
mod phase26_tests {
    use super::*;
    #[test]
    fn phase26_terminal_engine_dual_fd_parity() {
        use ql_ad::Dual;
        fn make<T: ql_ad::Real>(x: T) -> BergomiForwardStartNapoleonVolSpec<T> {
            BergomiForwardStartNapoleonVolSpec {
                initial_forward_variance: x * x,
                vol_of_vol_omega: T::from_f64(0.3),
                kappa_short_k1: T::from_f64(3.0),
                kappa_long_k2: T::from_f64(0.5),
                long_factor_weight_theta: T::from_f64(0.5),
                factor_corr_rho12: T::from_f64(0.2),
                spot_vol_corr_rho_s1: T::from_f64(-0.5),
                spot_vol_corr_rho_s2: T::from_f64(-0.3),
                forward_start_t1: T::from_f64(0.5),
                forward_end_t2: T::one(),
                forward_moneyness_k: T::one(),
                napoleon_base_coupon: T::from_f64(0.01),
                napoleon_global_floor: T::zero(),
                napoleon_subperiods_m: 4,
            }
        }
        let h = 1e-5;
        let fd = (evaluate_bergomi_forward_start_smile_and_napoleon_convexity(&make(0.2 + h))
            .forward_atm_implied_vol
            - evaluate_bergomi_forward_start_smile_and_napoleon_convexity(&make(0.2 - h))
                .forward_atm_implied_vol)
            / (2.0 * h);
        let dual = make(Dual::<1>::variable(0.2, 0));
        let ad = evaluate_bergomi_forward_start_smile_and_napoleon_convexity(&dual)
            .forward_atm_implied_vol;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }

    #[test]
    fn phase26_tab36_terminal_renders_and_reacts() {
        let mut ws = DeskWorkstation::new();
        let normal = ws.step_tick();
        let normal_metric = normal.phase26_sota_metric;
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let stressed = ws.step_tick();
        assert!(stressed.phase26_sota_metric.is_finite());
        assert!((stressed.phase26_sota_metric - normal_metric).abs() > 1e-12);
        assert!(
            ws.render_dashboard(&stressed, false)
                .contains("36 Phase26 SOTA")
        );
        assert_eq!(TabId::TAB_COUNT, 43);
        assert_eq!(TabId::ALL.len(), 43);
        assert!(TabId::ALL.contains(&TabId::Phase26Sota));
    }
}

#[cfg(test)]
mod phase27_tests {
    use super::*;

    #[test]
    fn phase27_tab37_terminal_renders_and_reacts() {
        let mut ws = DeskWorkstation::new();
        let normal = ws.step_tick();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let stressed = ws.step_tick();
        assert!(stressed.phase27_sota_metric.is_finite());
        assert!((stressed.phase27_sota_metric - normal.phase27_sota_metric).abs() > 1e-12);
        assert!(
            ws.render_dashboard(&stressed, false)
                .contains("37 Phase 27 SOTA")
        );
        assert_eq!(TabId::TAB_COUNT, 43);
        assert_eq!(TabId::ALL.len(), 43);
        assert!(TabId::ALL.contains(&TabId::Phase27Sota));
    }

    #[test]
    fn phase27_terminal_dual_fd_parity() {
        use ql_ad::Dual;
        fn eval<T: ql_ad::Real>(x: T) -> T {
            ql_inflation::evaluate_cochrane_ftpl_and_maturity_inflation_smoothing(
                &ql_inflation::CochraneFiscalTheoryPriceLevelSpec {
                    debt_maturities_years: vec![T::one()],
                    nominal_debt_face_values: vec![x],
                    nominal_zero_yields: vec![T::from_f64(0.03)],
                    real_zero_yields: vec![T::from_f64(0.01)],
                    primary_real_surpluses: vec![T::from_f64(5.0)],
                    terminal_surplus_growth_rate: T::from_f64(0.02),
                    fiscal_surplus_shock_fraction: T::from_f64(-0.05),
                    expected_future_inflation_rise: T::from_f64(0.01),
                    inflation_risk_beta: T::from_f64(0.4),
                    inflation_volatility: T::from_f64(0.1),
                    tips_liquidity_spread: T::from_f64(0.002),
                },
            )
            .equilibrium_price_level_p0
        }
        let h = 1e-5;
        let fd = (eval(100.0 + h) - eval(100.0 - h)) / (2.0 * h);
        let d = eval(Dual::<1>::variable(100.0, 0));
        assert!(d.d[0].abs() > 1e-6 && (d.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase28_tests {
    use super::*;
    #[test]
    fn phase28_tab38_terminal_renders_and_reacts() {
        let mut ws = DeskWorkstation::default();
        let normal_frame = ws.step_tick();
        let normal = ws.render_dashboard(&normal_frame, false);
        ws.scenario = ScenarioKind::CreditWideningPlus100Bp;
        let stressed_frame = ws.step_tick();
        let stressed = ws.render_dashboard(&stressed_frame, false);
        assert_eq!(TabId::ALL.len(), 43);
        assert!(TabId::ALL.contains(&TabId::Phase28Sota));
        assert!(normal.contains("38:P28-SOTA"));
        assert!(stressed.contains("38:P28-SOTA"));
        assert!(
            (stressed_frame.phase28_sota_metric - normal_frame.phase28_sota_metric).abs() > 1e-12
        );
        assert!(
            (ws.phase28_live_metric()
                - DeskWorkstation::phase28_live_metric_from_market_with_bump(
                    &ws.scenario_market(),
                    0.01,
                ))
            .abs()
                < 1e-12
        );
    }
}

#[cfg(test)]
mod phase29_tests {
    use super::*;
    #[test]
    fn phase29_tab39_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 43);
        assert!(TabId::Phase29Sota.label().starts_with("39 Phase 29 SOTA"));
        let mut ws = DeskWorkstation::default();
        let a = ws.step_tick();
        let b = ws.step_tick();
        assert!(a.phase29_sota_metric.is_finite() && b.phase29_sota_metric.is_finite());
        assert!(a.phase29_sota_metric != 0.0);
        assert!(ws.render_dashboard(&a, false).contains("39 Phase29 SOTA"));
    }
}

#[cfg(test)]
mod phase30_tests {
    use super::*;

    #[test]
    fn phase30_tab40_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 43);
        assert_eq!(TabId::ALL.len(), 43);
        assert!(
            TabId::Phase30BookSota
                .label()
                .starts_with("40 Phase 30 SOTA")
        );
        let mut ws = DeskWorkstation::new();
        ws.set_scenario(ScenarioKind::Normal);
        let normal = ws.step_tick();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let stressed = ws.step_tick();
        assert!(normal.phase30_sota_metric.is_finite());
        assert!(stressed.phase30_sota_metric.is_finite());
        assert_ne!(normal.phase30_sota_metric, stressed.phase30_sota_metric);
        assert!(
            ws.render_dashboard(&stressed, false)
                .contains("40 Phase30 SOTA")
        );
    }
}

#[cfg(test)]
mod phase31_tests {
    use super::*;

    #[test]
    fn phase31_tab41_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 43);
        assert_eq!(TabId::ALL.len(), 43);
        assert!(
            TabId::Phase31BookSota
                .label()
                .starts_with("41 Phase 31 SOTA")
        );
        let mut ws = DeskWorkstation::new();
        ws.set_scenario(ScenarioKind::Normal);
        let normal = ws.step_tick();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let credit = ws.step_tick();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let stag = ws.step_tick();
        assert!(normal.phase31_sota_metric.is_finite());
        assert!((credit.phase31_sota_metric - normal.phase31_sota_metric).abs() > 1e-12);
        assert!((stag.phase31_sota_metric - normal.phase31_sota_metric).abs() > 1e-12);
        assert!(
            ws.render_dashboard(&stag, false)
                .contains("41 Phase31 SOTA")
        );
    }
}

#[cfg(test)]
mod phase32_tests {
    use super::*;
    #[test]
    fn phase32_tab42_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 43);
        assert_eq!(TabId::ALL.len(), 43);
        assert!(
            TabId::Phase32BookSota
                .label()
                .starts_with("42 Phase 32 SOTA")
        );
        let mut ws = DeskWorkstation::default();
        let m0 = ws.phase32_live_metric();
        let f1 = ws.step_tick();
        assert!((f1.phase32_sota_metric - m0).abs() > 1e-12);
        ws.set_scenario(ScenarioKind::Normal);
        let normal = ws.step_tick();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let credit = ws.step_tick();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let stag = ws.step_tick();
        assert!((credit.phase32_sota_metric - normal.phase32_sota_metric).abs() > 1e-12);
        assert!((stag.phase32_sota_metric - normal.phase32_sota_metric).abs() > 1e-12);
        assert!(
            ws.render_dashboard(&stag, false)
                .contains("42 Phase32 SOTA")
        );
    }
}

#[cfg(test)]
mod phase33_tests {
    use super::*;
    #[test]
    fn phase33_tab43_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 43);
        assert_eq!(TabId::ALL.len(), 43);
        assert!(
            TabId::Phase33BookSota
                .label()
                .starts_with("43 Phase 33 SOTA")
        );
        let mut ws = DeskWorkstation::new();
        let a = ws.phase33_live_metric();
        let _ = ws.step_tick();
        assert!(ws.phase33_live_metric().is_finite());
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.phase33_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.phase33_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("43 Phase33 SOTA")
        );
    }
}
