import torch
import alphadeep.model.base as model_base
class Model(torch.nn.Module):
    def __init__(self,
        dropout=0.2
    ):
        super().__init__()
        self.aa_embedding_size = 27

        self.dropout = torch.nn.Dropout(dropout)

        hidden = 256
        self.rt_encoder = model_base.Input_AA_CNN_LSTM_Encoder(
            hidden
        )

        self.rt_decoder = model_base.LinearDecoder(
            hidden,
            1
        )

    def forward(self,
        aa_indices,
        mod_x,
    ):
        x = self.rt_encoder(aa_indices, mod_x)
        x = self.dropout(x)

        return self.rt_decoder(x).squeeze(1)
