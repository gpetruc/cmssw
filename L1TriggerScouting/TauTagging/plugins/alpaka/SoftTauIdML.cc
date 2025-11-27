#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftTauDeviceTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/AssociationMapDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/ClustersDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "L1TriggerScouting/TauTagging/plugins/alpaka/TransformKernel.h"
#include "PhysicsTools/PyTorchAlpaka/interface/QueueGuard.h"
#include "PhysicsTools/PyTorchAlpaka/interface/TensorRegistry.h"
#include "PhysicsTools/PyTorchAlpaka/interface/alpaka/AlpakaModel.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  class SoftTauIdML : public stream::EDProducer<> {
  public:
    SoftTauIdML(const edm::ParameterSet &params)
        : EDProducer<>(params),
          pf_token_(consumes(params.getParameter<edm::InputTag>("pf"))),
          association_map_token_{consumes(params.getParameter<edm::InputTag>("clusters"))},
          soft_tau_token_{produces()},
          model_(params.getParameter<edm::FileInPath>("model").fullPath()),
          do_inference_{params.getParameter<bool>("do_inference")},
          max_batch_size_{params.getParameter<uint32_t>("maxBatchSize")} {}

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::FileInPath>("model");
      desc.add<edm::InputTag>("pf");
      desc.add<edm::InputTag>("clusters");
      desc.add<bool>("do_inference");
      desc.add<uint32_t>("maxBatchSize", std::numeric_limits<uint32_t>::max());
      descriptions.addWithDefaultLabel(desc);
    }

    void produce(device::Event &event, const device::EventSetup &event_setup) override {
      // in/out collections
      const auto &pf = event.get(pf_token_);
      const auto &association_map = event.get(association_map_token_);
      SoftTauInputDeviceTensor input_tensor = kernels::transform(event.queue(), pf, association_map);

      // prepare output and initialize to zero
      const auto job_size = input_tensor.view().metadata().size();
      auto output_tensor = SoftTauOutputDeviceTensor(job_size, event.queue());
      output_tensor.zeroInitialise(event.queue());

      if (do_inference_) {
        const auto batch_size = std::min<uint32_t>(job_size, max_batch_size_);
  
        // records
        auto input_records = input_tensor.view().records();
        auto output_records = output_tensor.view().records();
        // input tensor definition
        cms::torch::alpakatools::TensorRegistry<Device> inputs(batch_size);
        inputs.register_tensor<SoftTauInputTensorSoA>("jet_features", input_records.features());
        inputs.register_tensor<SoftTauInputTensorSoA>("padding_mask", input_records.pad_mask());
        // output tensor definition
        cms::torch::alpakatools::TensorRegistry<Device> outputs(batch_size);
        outputs.register_tensor<SoftTauOutputTensorSoA>("cls", output_records.cls());
        outputs.register_tensor<SoftTauOutputTensorSoA>("reg_vz", output_records.vz());
        outputs.register_tensor<SoftTauOutputTensorSoA>("reg_pt", output_records.pt());
        outputs.register_tensor<SoftTauOutputTensorSoA>("charge", output_records.charge());
  
        // inference, queue guard restore stream when goes out of scope
        {
          cms::torch::alpakatools::QueueGuard<Queue> guard(event.queue());
          model_.to(event.queue());
          model_.forward(event.queue(), inputs, outputs);
        }
      }
      
      // put device-side product into event
      event.emplace(soft_tau_token_, std::move(output_tensor));
    }

  private:
    // event query tokens
    const device::EDGetToken<PFCandidateDeviceCollection> pf_token_;
    // clustering output
    const device::EDGetToken<AssociationMapDevice> association_map_token_;
    // put ml output into event
    const device::EDPutToken<SoftTauOutputDeviceTensor> soft_tau_token_;
    // model
    torch::AlpakaModel model_;
    // do inference or not
    const bool do_inference_;
    // scouting switch
    const uint32_t max_batch_size_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::torchtest

DEFINE_FWK_ALPAKA_MODULE(l1sc::SoftTauIdML);