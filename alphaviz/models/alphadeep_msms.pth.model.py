import torch
import alphadeep.model.base as model_base
class Model(torch.nn.Module):
    def __init__(self,
        num_frag_types,
        num_modloss_types=0,
        mask_modloss=True,
        dropout=0.2,
    ):
        super().__init__()

        self.dropout = torch.nn.Dropout(dropout)

        BiRNN = True
        hidden=512
        hidden_rnn_layer=2

        self.input_nn = model_base.InputAALSTM_cat_Meta(hidden)

        self.hidden_nn = model_base.SeqLSTM(
            hidden, hidden, rnn_layer=hidden_rnn_layer,
            bidirectional=BiRNN
        )

        self.output_nn = model_base.OutputLSTM_cat_Meta(
            hidden,
            num_frag_types,
        )

        self._num_modloss_types = num_modloss_types
        if num_modloss_types:
            # for transfer learning of modloss frags
            self.modloss_nn = torch.nn.ModuleList([
                model_base.SeqLSTM(
                    hidden, hidden,
                    rnn_layer=1, bidirectional=BiRNN
                ),
                model_base.SeqLSTM(
                    hidden, num_modloss_types,
                    rnn_layer=1, bidirectional=False
                ),
            ])
            self._mask_modloss = mask_modloss
            self._non_modloss = num_frag_types-num_modloss_types
        else:
            self._mask_modloss = True

    def forward(self,
        aa_indices,
        mod_x,
        charges:torch.Tensor,
        NCEs:torch.Tensor,
        instrument_indices,
    ):

        in_x = self.input_nn(
            aa_indices, mod_x,
            charges, NCEs, instrument_indices
        )
        in_x = self.dropout(in_x)

        hidden_x = self.hidden_nn(in_x)
        hidden_x = self.dropout(hidden_x)

        out_x = self.output_nn(
            hidden_x,
            charges, NCEs, instrument_indices
        )

        # modloss is mainly only for Phospho@S/T
        if not self._mask_modloss:
            modloss_x = self.modloss_nn[0](
                in_x
            )*in_x + hidden_x
            modloss_x = self.modloss_nn[-1](
                modloss_x
            )
            out_x = torch.cat((
                out_x[:,:,:self._non_modloss],
                out_x[:,:,self._non_modloss:]+modloss_x
            ),2)

        return out_x[:,3:,:]
